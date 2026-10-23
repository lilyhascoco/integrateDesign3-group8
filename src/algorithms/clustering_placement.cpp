#include "algorithms/clustering_placement.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace plc {
namespace {
struct Cluster {
    std::vector<CellId> cells;
    std::map<CellId, double> edges; // ordered traversal makes tie handling reproducible
    double area = 0, x = 0, y = 0;
    unsigned version = 0;
    bool active = false;
};
struct Edge {
    double weight;
    CellId a, b;
    unsigned va, vb;
    bool operator<(const Edge& e) const {
        if (weight != e.weight) return weight < e.weight;
        return std::tie(a, b) > std::tie(e.a, e.b);
    }
};
}
void ClusteringPlacement::run(const PlacementDB& db, PlacementState* state) {
    stats_ = {};
    if (!state) throw std::invalid_argument("Null placement state");
    state->resize(db.num_cells());
    const auto& core = db.design.core;
    const double w = core.ux - core.lx, h = core.uy - core.ly;
    if (!(w > 0 && h > 0)) throw std::invalid_argument("Invalid core region");
    std::vector<Cluster> c(db.num_cells());
    std::vector<CellId> owner(db.num_cells());
    for (const auto& cell : db.design.cells) {
        state->set_pos(cell.id, cell.init_x, cell.init_y);
        state->set_orient(cell.id, cell.orient);
        owner[cell.id] = cell.id;
        if (!cell.is_movable()) continue;
        ++stats_.movable_cells;
        if (cell.width > w || cell.height > h || cell.width < 0 || cell.height < 0) {
            ++stats_.skipped_cells;
            continue;
        }
        auto& a = c[cell.id];
        a.active = true;
        a.cells.push_back(cell.id);
        a.area = cell.area();
        a.x = std::clamp(cell.init_x + cell.width/2, core.lx, core.ux);
        a.y = std::clamp(cell.init_y + cell.height/2, core.ly, core.uy);
        ++stats_.clusters;
    }
    auto connect = [&](CellId a, CellId b, double weight) {
        if (a != b && c[a].active && c[b].active) {
            c[a].edges[b] += weight;
            c[b].edges[a] += weight;
        }
    };
    for (const auto& net : db.design.nets) {
        if (!(net.weight > 0)) continue;
        std::vector<CellId> ids;
        for (auto p : net.pin_ids) ids.push_back(db.pin(p).cell_id);
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        if (ids.size() < 2) continue;
        const double weight = net.weight / (ids.size()-1);
        if (ids.size() <= std::max<std::size_t>(2, config_.clique_limit)) {
            for (std::size_t i=0;i<ids.size();++i)
                for (std::size_t j=i+1;j<ids.size();++j) connect(ids[i],ids[j],weight);
        } else {
            // Deterministic ring sparsification: retain all endpoints, O(degree) edges.
            ++stats_.approximated_nets;
            for (std::size_t i=0;i<ids.size();++i)
                connect(ids[i],ids[(i+1)%ids.size()],net.weight/2);
        }
    }
    std::priority_queue<Edge> heap;
    const auto cap = std::max<std::size_t>(1, config_.max_cluster_cells);
    auto push = [&](CellId a, CellId b, double weight) {
        if (a > b) std::swap(a,b);
        if (c[a].cells.size()+c[b].cells.size() <= cap)
            heap.push({weight,a,b,c[a].version,c[b].version});
    };
    for (CellId i=0;i<c.size();++i)
        for (const auto& e : c[i].edges) if (i<e.first) push(i,e.first,e.second);
    const auto target = config_.target_clusters ? config_.target_clusters :
        std::max<std::size_t>(1,(stats_.clusters+31)/32);
    while (stats_.clusters > target && !heap.empty()) {
        const auto e=heap.top(); heap.pop();
        auto& a=c[e.a]; auto& b=c[e.b];
        if (!a.active || !b.active || a.version!=e.va || b.version!=e.vb) continue;
        if (a.cells.size()+b.cells.size()>cap) continue;
        double aw=a.area, bw=b.area;
        if (aw+bw==0) { aw=a.cells.size(); bw=b.cells.size(); }
        a.x=(a.x*aw+b.x*bw)/(aw+bw); a.y=(a.y*aw+b.y*bw)/(aw+bw);
        a.area+=b.area;
        a.cells.insert(a.cells.end(),b.cells.begin(),b.cells.end());
        for (auto id:b.cells) owner[id]=e.a;
        a.edges.erase(e.b);
        for (const auto& [neighbor,weight]:b.edges) {
            if (neighbor==e.a) continue;
            c[neighbor].edges.erase(e.b);
            a.edges[neighbor]+=weight;
            c[neighbor].edges[e.a]=a.edges[neighbor];
        }
        b.active=false; b.edges.clear(); b.cells.clear(); ++b.version; ++a.version;
        for (const auto& [neighbor,weight]:a.edges) push(e.a,neighbor,weight);
        --stats_.clusters; ++stats_.merges;
    }
    std::vector<CellId> active;
    for (CellId i=0;i<c.size();++i) if(c[i].active) active.push_back(i);
    if(active.empty()) return;
    // Spread coarse clusters before connectivity relaxation; the tether prevents collapse.
    const auto cols=std::max<std::size_t>(1,std::ceil(std::sqrt(active.size()*w/h)));
    const auto rows=(active.size()+cols-1)/cols;
    std::vector<double> anchor_x(c.size()),anchor_y(c.size());
    for(std::size_t k=0;k<active.size();++k) {
        auto id=active[k];
        c[id].x=anchor_x[id]=core.lx+(k%cols+0.5)*w/cols;
        c[id].y=anchor_y[id]=core.ly+(k/cols+0.5)*h/rows;
    }
    // Fixed terminals contribute anchors through incident nets, never enter movable clusters.
    std::vector<double> fx(c.size()),fy(c.size()),fw(c.size());
    for(const auto& net:db.design.nets) {
        double x=0,y=0; std::size_t count=0;
        for(auto p:net.pin_ids) {
            const auto& pin=db.pin(p); const auto& cell=db.cell(pin.cell_id);
            if(cell.is_fixed()) {x+=cell.init_x+cell.width/2+pin.offset_x;
                y+=cell.init_y+cell.height/2+pin.offset_y; ++count;}
        }
        if(!count || !(net.weight>0)) continue;
        for(auto p:net.pin_ids) {
            auto id=owner[db.pin(p).cell_id];
            if(c[id].active) {fx[id]+=net.weight*x/count; fy[id]+=net.weight*y/count; fw[id]+=net.weight;}
        }
    }
    std::vector<double> nx(c.size()),ny(c.size());
    for(std::size_t iter=0;iter<config_.coarse_iterations;++iter) {
        for(auto id:active) {
            double sum=fw[id],x=fx[id],y=fy[id];
            for(const auto& [j,weight]:c[id].edges) {sum+=weight; x+=weight*c[j].x; y+=weight*c[j].y;}
            const double tether=std::max(1.0,sum*0.1);
            nx[id]=std::clamp((x+tether*anchor_x[id])/(sum+tether),core.lx,core.ux);
            ny[id]=std::clamp((y+tether*anchor_y[id])/(sum+tether),core.ly,core.uy);
        }
        for(auto id:active) {c[id].x=nx[id];c[id].y=ny[id];}
    }
    for(auto id:active) {
        auto& a=c[id]; std::sort(a.cells.begin(),a.cells.end());
        const double boxw=std::min(w,std::sqrt(std::max(1.0,a.area)*w/h));
        const double boxh=std::min(h,std::max(1.0,a.area)/boxw);
        const auto nc=std::max<std::size_t>(1,std::ceil(std::sqrt(a.cells.size()*boxw/boxh)));
        const auto nr=(a.cells.size()+nc-1)/nc;
        for(std::size_t k=0;k<a.cells.size();++k) {
            const auto& cell=db.cell(a.cells[k]);
            double x=a.x-boxw/2+(k%nc+0.5)*boxw/nc-cell.width/2;
            double y=a.y-boxh/2+(k/nc+0.5)*boxh/nr-cell.height/2;
            state->set_pos(cell.id,std::clamp(x,core.lx,core.ux-cell.width),
                std::clamp(y,core.ly,core.uy-cell.height));
            ++stats_.placed_cells;
        }
    }
}
} // namespace plc
