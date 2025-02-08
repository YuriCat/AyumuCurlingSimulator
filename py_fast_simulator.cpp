
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "fast_simulator.hpp"


using namespace std;
using namespace FastSimulator;

namespace py = pybind11;


vector<tuple<double, double>> py_simulate(vector<tuple<double, double>> stones, int index, tuple<double, double, int> shot, bool freeguard, bool rink_only) {
    Sheet sheet;
    sheet.clear();
    int count = 0;
    for (auto stone : stones) {
        auto [x, y] = stone;
        if (y > 0) sheet.set_stone(count, Position(x, y)); // check valid stone by y > 0
        count += 1;
    }
    auto [vx, vy, spin] = shot;
    bool foul = simulate(sheet, index, Shot(vx, vy, spin), freeguard, rink_only);
    vector<tuple<double, double>> results(16, make_tuple(0, 0));
    if (foul) {
        int count = 0;
        for (auto stone : stones) {
            results[count] = stone;
            count += 1;
        }
    } else {
        for (int i = 0; i < 16; i++) {
            if (sheet.stone_bits & (1U << i)) results[i] = make_tuple(sheet.stone[i].x, sheet.stone[i].y);
        }
    }
    return results;
}

tuple<double, double> py_shot2dest(tuple<double, double, int> shot) {
    auto [vx, vy, spin] = shot;
    Position p = shot2dest(Shot(vx, vy, spin));
    return make_tuple(p.x, p.y);
}

tuple<double, double, int> py_dest2shot(tuple<double, double> p, int spin) {
    auto [x, y] = p;
    Shot shot = dest2shot(Position(x, y), spin);
    return make_tuple(shot.vx, shot.vy, spin);
}

tuple<double, double, int> py_passpoint2shot(tuple<double, double> p, double v, int spin) {
    auto [x, y] = p;
    Shot shot = passpoint2shot(Position(x, y), v, spin);
    return make_tuple(shot.vx, shot.vy, spin);
}

tuple<double, double, int> py_passpointgo2shot(tuple<double, double> p, double r, int spin) {
    auto [x, y] = p;
    Shot shot = passpointgo2shot(Position(x, y), r, spin);
    return make_tuple(shot.vx, shot.vy, spin);
}


PYBIND11_MODULE(fast_simulator, m)
{
    m.doc() = "digital curling fast simulator";

    m.def("simulate", &py_simulate, "simulation", py::arg("stones"), py::arg("index"), py::arg("shot"), py::arg("freeguard") = false, py::arg("rink_only") = true);

    m.def("shot2dest", &py_shot2dest, "compute destination from shot");
    m.def("dest2shot", &py_dest2shot, "compute shot from destination");
    m.def("passpoint2shot", &py_passpoint2shot, "compute shot from passing point and velocity");
    m.def("passpointgo2shot", &py_passpointgo2shot, "compute shot from passing point and distance after that point");
};
