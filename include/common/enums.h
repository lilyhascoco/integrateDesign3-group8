#pragma once

namespace plc {

enum class CellType {
    kStdCell,
    kMacro,
    kTerminal,
    kTerminalNI,
    kUnknown
};

enum class PlaceStatus {
    kMovable,
    kFixed,
    kFixedNI,
    kUnplaced
};

enum class Orient {
    kN,
    kS,
    kE,
    kW,
    kFN,
    kFS,
    kFE,
    kFW,
    kUnknown
};

enum class PinDirection {
    kInput,
    kOutput,
    kInout,
    kUnknown
};

}  // namespace plc
