
#include <iostream>
#include "fast_simulator.hpp"

using namespace std;
using namespace FastSimulator;

int main() {
    // 0: right spin  1: left spin
    Position pos;
    Shot shot;
    Sheet sheet;

    pos = shot2dest(Shot(-0.132, 2.3995, 0));
    cout << pos << endl;

    shot = dest2shot(Position(0, 38.405), 0);
    cout << shot << endl;

    shot = passpoint2shot(Position(0, 38.405), 3.5, 0);
    cout << shot << endl;

    shot = passpointgo2shot(Position(0, 38.405), 2.0, 1);
    cout << shot << endl;

    sheet.clear();
    simulate(sheet, 15, Shot(-0.132, 2.3995, 0));
    cout << sheet << endl;

    // hitting simulation
    sheet.clear();
    sheet.set_stone(0, Position(0, 38.405));
    shot = passpoint2shot(sheet.stone[0], 3.5, 1);
    simulate(sheet, 1, shot, false);
    cout << sheet << endl;

    sheet.clear();
    sheet.set_stone(0, Position(0, 38.405));
    simulate(sheet, 1, shot, false, false);
    cout << sheet << endl;

    // pushing simulation
    sheet.clear();
    sheet.set_stone(0, Position(0, 36));
    shot = passpointgo2shot(sheet.stone[0], 2.4, 0);
    simulate(sheet, 1, shot);
    cout << sheet << endl;

    // freeguard foul
    sheet.clear();
    sheet.set_stone(0, Position(0, 36));
    shot = passpoint2shot(sheet.stone[0], 3.5, 1);
    simulate(sheet, 1, shot);
    cout << sheet << endl;

    sheet.clear();
    sheet.set_stone(0, Position(0, 36));
    simulate(sheet, 1, shot, true);
    cout << sheet << endl;

    return 0;
}
