
import fast_simulator as fs


# 0: right spin  1: left spin

pos = fs.shot2dest((-0.132, 2.3995, 0))
print(pos)

shot = fs.dest2shot((0, 38.405), 0)
print(shot)

shot = fs.passpoint2shot((0, 38.405), 3.5, 0)
print(shot)

shot = fs.passpointgo2shot((0, 38.405), 2.0, 1)
print(shot)

stones_solo = fs.simulate([], 15, (-0.132, 2.3995, 0))
print(stones_solo)

# hitting simulation
stones = [(0, 38.405)]
shot = fs.passpoint2shot(stones[0], 3.5, 1)
stones_rink = fs.simulate(stones, 1, shot, freeguard=False)
stones_plane = fs.simulate(stones, 1, shot, freeguard=False, rink_only=False)

print(stones_rink[:2])
print(stones_plane[:2])

# pushing simulation
stones = [(0, 36)]
shot = fs.passpointgo2shot(stones[0], 2.4, 0)
stones_ = fs.simulate(stones, 1, shot)

print(stones_[:2])

# freeguard foul
stones = [(0, 36)]
shot = fs.passpoint2shot(stones[0], 3.5, 1)
stones_ = fs.simulate(stones, 1, shot)
stones_fg = fs.simulate(stones, 1, shot, freeguard=True)

print(stones_[:2])
print(stones_fg[:2])
