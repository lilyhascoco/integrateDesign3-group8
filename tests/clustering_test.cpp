#include "algorithms/quadratic_analytical_placement.h"
#include "algorithms/local_sa_placer.h"
#include "algorithms/clustering_placement.h"
#include "metrics/hpwl.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
void check(bool v) {if(!v) throw std::runtime_error("Test assertion failed");}
int main() {
    plc::PlacementDB db; db.design.core={100,200,300,400};
    for(unsigned i=0;i<5;++i) {
        plc::Cell c; c.id=i;c.width=10+i*2;c.height=10;c.init_x=120;c.init_y=220;
        c.place_status=i==4?plc::PlaceStatus::kFixed:plc::PlaceStatus::kMovable;
        c.orient=i==4?plc::Orient::kFN:plc::Orient::kN;
        db.design.cells.push_back(c);
    }
    auto edge=[&](unsigned a,unsigned b,double weight) {
        plc::Net n;n.id=db.design.nets.size();n.weight=weight;
        for(auto id:{a,b}) {plc::Pin p;p.cell_id=id;p.id=db.design.pins.size();p.net_id=n.id;
            n.pin_ids.push_back(p.id);db.design.pins.push_back(p);}
        db.design.nets.push_back(n);
    };
    edge(0,1,10);edge(1,2,1);edge(2,3,9);edge(3,4,1);
    plc::ClusteringConfig cfg;cfg.target_clusters=2;cfg.max_cluster_cells=2;
    plc::ClusteringPlacement placer(cfg);plc::PlacementState a,b;
    placer.run(db,&a);check(placer.stats().clusters==2);check(placer.stats().merges==2);
    check(placer.stats().placed_cells==4);
    placer.run(db,&b);
    for(auto& c:db.design.cells) {
        check(a.x(c.id)==b.x(c.id)&&a.y(c.id)==b.y(c.id));
        if(c.is_fixed()) check(a.x(c.id)==c.init_x && a.y(c.id)==c.init_y);
        else check(a.x(c.id)>=100&&a.y(c.id)>=200&&a.x(c.id)+c.width<=300&&a.y(c.id)+c.height<=400);
    }
    auto distance=[&](unsigned i,unsigned j) {return std::hypot(a.x(i)-a.x(j),a.y(i)-a.y(j));};
    check(distance(0,1)<distance(0,2));check(distance(2,3)<distance(0,2));
    plc::PlacementState quadratic_state;
    plc::QuadraticAnalyticalPlacement quadratic;
    quadratic.run(db, &quadratic_state);
    check(quadratic_state.x(4)==db.cell(4).init_x);
    check(quadratic_state.y(4)==db.cell(4).init_y);
    check(quadratic_state.orient(4)==db.cell(4).orient);
    plc::LocalSAConfig sa_config;sa_config.bins_x=1;sa_config.bins_y=1;sa_config.max_iterations=20;
    const double before_sa=plc::HpwlEvaluator{}.compute_total_hpwl(db, a);
    plc::LocalSAPlacer sa(sa_config,42);sa.run(db,&a);
    check(std::abs(sa.stats().initial_hpwl-before_sa)<1e-8);
    check(std::abs(sa.stats().final_hpwl-plc::HpwlEvaluator{}.compute_total_hpwl(db,a))<1e-8);
    // Known center-relative pin coordinates: centers (5,5) and (22,5), weight 2 -> 34.
    plc::PlacementDB hd;hd.design.cells.resize(2);hd.design.cells[0].width=10;hd.design.cells[0].height=10;
    hd.design.cells[1].width=20;hd.design.cells[1].height=10;
    hd.design.pins.resize(2);hd.design.pins[1].cell_id=1;
    plc::Net n;n.weight=2;n.pin_ids={0,1};hd.design.nets.push_back(n);
    plc::PlacementState hs;hs.resize(2);hs.set_pos(0,0,0);hs.set_pos(1,12,0);
    check(plc::HpwlEvaluator{}.compute_total_hpwl(hd,hs)==34);
    cfg.clique_limit=2;plc::Net high;high.pin_ids={0,1,2,3};db.design.nets.push_back(high);
    plc::ClusteringPlacement approx(cfg);approx.run(db,&a);check(approx.stats().approximated_nets==1);
    db.design.nets.clear();placer.run(db,&a);check(placer.stats().clusters==4);
    db.design.cells[0].width=201;placer.run(db,&a);check(placer.stats().skipped_cells==1);
    db.design.cells.clear();placer.run(db,&a);check(placer.stats().clusters==0);
    std::cout<<"Clustering and HPWL tests passed\n";
}
