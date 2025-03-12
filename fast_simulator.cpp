
// Ayumu Curling Simulator
// Katsuki Ohto

#include <array>
#include <tuple>
#include <vector>
#include "fast_simulator.hpp"


namespace FastSimulator {
    float traj_table[__TABLE + 1][7]; // 0 <= t <= 100
    float v2rtt_table[__TABLE + 1][4]; // 0 <= v <= 5
    float r2vt_table[__TABLE + 1][2]; // 0 <= r <= 196

    Initializer __initializer;

    // solo simulation
    PositionTime simulate_solo(Shot shot, fpn_t time_step) {
        fpn_t x = X_THROW;
        fpn_t y = Y_THROW;
        fpn_t vx = shot.vx;
        fpn_t vy = shot.vy;
        fpn_t w = spin2wthrow(shot.spin);
        int t = 0;
        friction_step(vx, vy, w, &vx, &vy, &w, time_step / 2);
        while (1) {
            if (vx == 0 && vy == 0) break;
            x += vx * time_step;
            y += vy * time_step;
            friction_step(vx, vy, w, &vx, &vy, &w, time_step);
            t += 1;
        }
        return PositionTime(x, y, t * time_step);
    }

    // functions for fast simulation
    template <typename T>
    inline T clip(T v, T min, T max) { return std::min(max, std::max(min, v)); }

    void collision(Stone *s0, Stone *s1, bool stopped) {
        Angle coltheta = xy2ang(s1->x - s0->x, s1->y - s0->y);

        fpn_t v0n = s0->v * (s0->th - coltheta).c;
        fpn_t v0t = s0->v * (s0->th - coltheta).s;
        fpn_t v1n = stopped ? 0 : s1->v * (s1->th - coltheta).c;
        fpn_t v1t = stopped ? 0 : s1->v * (s1->th - coltheta).s;
        fpn_t w0 = s0->w;
        fpn_t w1 = stopped ? 0 : s1->w;

        fpn_t vn = v0n - v1n;
        fpn_t vt = v0t - v1t;

        fpn_t ni_m = vn;
        if (vn <= RESTITUTION_THRESHOLD) ni_m /= 2;
        fpn_t max_friction = FRICTION_STONES * ni_m;

        fpn_t old_tangent_lambda = fpn_t(1 / 6.0) * (vt - (w0 + w1) * R_STONE_RADIUS);
        fpn_t ti_m = clip(old_tangent_lambda, -max_friction, +max_friction);

        fpn_t dw = ti_m * fpn_t(2.0 / R_STONE_RADIUS);

        // update
        s0->v = std::min(V_SIMULATION_MAX, xy2r(v0n - ni_m, v0t - ti_m));
        s0->th = xy2ang(v0t - ti_m, v0n - ni_m) + coltheta;
        s0->w += dw;

        s1->v = std::min(V_SIMULATION_MAX, xy2r(v1n + ni_m, v1t + ti_m));
        s1->th = xy2ang(v1t + ti_m, v1n + ni_m) + coltheta;
        if (stopped) s1->w = dw; else s1->w += dw;
    }

    bool has_collision_chance_aw_as(const Stone& aw, const Stone& as) {
        // determine the possibility of collision between a moving stone and a stopped stone
        // 動く石と静止している石の衝突可能性判定
        fpn_t dx = as.x - aw.x;
        fpn_t dy = as.y - aw.y;
        if (xy2r2(dx, dy) > r2r2(aw.gr + 2 * R_STONE_RADIUS)) return false; // R judge
        Angle l, r;
        if (!is_left_spin(w2spin(aw.w))) { l = aw.th;  r = aw.gth; }
        else                             { l = aw.gth; r = aw.th;  }
        if (dx * l.c - dy * l.s < -2 * R_STONE_RADIUS) return false; // side judge
        if (dx * r.c - dy * r.s > +2 * R_STONE_RADIUS) return false; // side sudge
        if (dx * l.s + dy * l.c < 0) return false; // back judge
        if (dx * r.s + dy * r.c < 0) return false; // back judge
        return true;
    }

    fpn_t min_contact_mr_aw_as(const Stone& aw, const Stone& as) {
        // lower bound on distance to collision between moving stone and stopped stone
        // 静止している石に衝突するまでの距離の下界
        Angle l, r;
        if (!is_left_spin(w2spin(aw.w))) { l = aw.th;  r = aw.gth; }
        else                             { l = aw.gth; r = aw.th;  }
        fpn_t mr;
        fpn_t dx = as.x - aw.x;
        fpn_t dy = as.y - aw.y;
        fpn_t r2 = distance2(aw, as);
        if (dx * l.c - dy * l.s < 0) {
            fpn_t rk = dx * l.s + dy * l.c;
            mr = rk - std::sqrt(rk * rk - r2 + 4 * R_STONE_RADIUS * R_STONE_RADIUS);
        } else if (dx * r.c - dy * r.s > 0) {
            fpn_t rk = dx * r.s + dy * r.c;
            mr = rk - std::sqrt(rk * rk - r2 + 4 * R_STONE_RADIUS * R_STONE_RADIUS);
        } else mr = stone_gap(aw, as);
        return mr;
    }

    fpn_t min_contact_time_2aw(const Stone& aw0, const Stone& aw1, fpn_t basetime) {
        // lower bound on collision time between moving stones
        // 動いている石同士の衝突時刻の下界
        if (distance2(aw0, aw1) > r2r2(aw0.gr + aw1.gr + 2 * R_STONE_RADIUS)) return 100000;
        fpn_t r = stone_gap(aw0, aw1);
        fpn_t t = r / (aw0.v + aw1.v);
        if (t >= basetime) return basetime;

        Angle coltheta0 = xy2ang(aw1.x - aw0.x, aw1.y - aw0.y);
        Angle coltheta1 = coltheta0 - Angle(0, -1);

        Angle l0, r0, l1, r1;
        if (!is_left_spin(w2spin(aw0.w))) { l0 = aw0.th;  r0 = aw0.gth; }
        else                              { l0 = aw0.gth; r0 = aw0.th;  }
        if (!is_left_spin(w2spin(aw1.w))) { l1 = aw1.th;  r1 = aw1.gth; }
        else                              { l1 = aw1.gth; r1 = aw1.th;  }

        fpn_t c0, c1;
        if      (coltheta0 < l0) c0 = std::max(fpn_t(0), (coltheta0 - l0).c);
        else if (coltheta0 > r0) c0 = std::max(fpn_t(0), (coltheta0 - r0).c);
        else c0 = 1;
        if      (coltheta1 < l1) c1 = std::max(fpn_t(0), (coltheta1 - l1).c);
        else if (coltheta1 > r1) c1 = std::max(fpn_t(0), (coltheta1 - r1).c);
        else c1 = 1;

        fpn_t nr0 = aw0.gr * c0;
        fpn_t nr1 = aw1.gr * c1;
        if (nr0 + nr1 < r) return 100000;

        fpn_t nv0 = aw0.v * c0;
        fpn_t nv1 = aw1.v * c1;
        return nv0 + nv1 > 0 ? r / (nv0 + nv1) : 100000;
    }

    void Stone::set_goal_info() {
        if (v == 0) {
            gr = gt = wp = 0;
            gth = oth = Angle(0, 0);
            gx = x; gy = y;
            return;
        }
        gr = __v2r(v);
        Angle gdth = __v2th(v);
        gt = __v2t(v);
        oth = __t2gth(gt);
        wp = __t2wp(gt);
        if (is_left_spin(w2spin(w))) {
            gdth = -gdth;
            oth = -oth;
        } else {
            wp = -wp;
        }
        gth = th + gdth;
        fpn_t dx, dy;
        rotate(gr, fpn_t(0), gth, &dy, &dx);
        gx = x + dx;
        gy = y + dy;
    }

    void Stone::step_to_stop() {
        x = gx; y = gy;
        v = w = 0;
    }

    void Stone::step_by_next_time(fpn_t nt) {
        fpn_t ngt = gt - nt;
        fpn_t ngr = __t2r(ngt);
        Angle noth = __t2gth(ngt);
        fpn_t nv = __t2v(ngt);
        Angle novth = __t2vth(ngt);
        fpn_t nwp = __t2wp(ngt);
        if (std::signbit(oth.s)) {
            noth = -noth;
            novth = -novth;
            w = std::max(fpn_t(0), w - nwp + wp);
        } else {
            nwp = -nwp;
            w = std::min(-fpn_t(0), w - nwp + wp);
        }
        Angle ngth = gth + oth - noth;
        fpn_t dx, dy;
        rotate(-ngr, fpn_t(0), ngth, &dy, &dx);
        x = gx + dx;
        y = gy + dy;
        v = nv;
        th = gth + oth - novth;
        // step goal information
        gr = ngr;
        gth = ngth;
        gt = ngt;
        oth = noth;
        wp = nwp;
    }

    BitSet32 resolve_contacts_(Sheet& sheet, std::tuple<int, int, fpn_t, fpn_t> *contact, int num_contacts) {
        // resolve contacts until there is no heading contacts
        // 向き合っているコンタクトが無くなるまで繰り返し衝突を処理
        BitSet32 updated(0);
        while (1) {
            // find the earliest collision
            // 最も速く衝突するコンタクトを探す
            int min_index = -1;
            fpn_t min_time = 100000;
            for (int i = 0; i < num_contacts; i++) {
                auto [i0, i1, d, sp] = contact[i];
                if (sp > 0) {
                    fpn_t t = d / std::max(fpn_t(1e-16), sp);
                    if (t < min_time) {
                        min_index = i;
                        min_time = t;
                    }
                }
            }
            if (min_index < 0) break;

            // collision 衝突
            auto [index0, index1, _, __] = contact[min_index];
            collision(&sheet.stone[index0], &sheet.stone[index1]);
            updated.set(index0);
            updated.set(index1);

            // update distance and relative velocity 距離と相対速度の更新
            if (num_contacts == 1) break;
            for (int i = 0; i < num_contacts; i++) {
                auto [i0, i1, d, sp] = contact[i];
                if (i == min_index) {
                    contact[i] = std::make_tuple(i0, i1, fpn_t(0), -sp);
                } else if (i0 == index0 || i0 == index1 || i1 == index0 || i1 == index1) {
                    contact[i] = std::make_tuple(i0, i1, d - min_time * sp, heading_speed(sheet.stone[i0], sheet.stone[i1]));
                } else if (sp != 0) {
                    contact[i] = std::make_tuple(i0, i1, d - min_time * sp, sp);
                }
            }
        }

        sheet.moving_stone_bits |= updated;
        return updated;
    }

    BitSet32 resolve_contacts(Sheet& sheet) {
        // check contacts
        // コンタクトの状況を調べる
        std::tuple<int, int, fpn_t, fpn_t> contact[32];
        int num_contacts = 0;
        for (int i = 0; i < 16; i++) {
            if (!(sheet.stone_bits & (1 << i))) continue;
            for (int j = 0; j < i; j++) {
                if (!(sheet.stone_bits & (1 << j))) continue;
                const Stone& s0 = sheet.stone[i];
                const Stone& s1 = sheet.stone[j];
                if (has_contact(s0, s1)) contact[num_contacts++] =
                    std::make_tuple(i, j, distance(s0, s1), heading_speed(s0, s1));
            }
        }
        BitSet32 updated = resolve_contacts_(sheet, contact, num_contacts);

        sheet.moving_stone_bits |= updated;
        return updated;
    }

    // fast simulation
    bool simulate(Sheet& sheet, bool rink_only) {
        BitSet32 awake = sheet.moving_stone_bits; // moving stone bits
        BitSet32 asleep = sheet.stone_bits & ~awake; // stopped stone bits
        BitSet32 collidable[16]; //collidable stone bits from each moving stones それぞれの動石から衝突可能な石のビット集合
        std::tuple<int, int, fpn_t, fpn_t> contact[32];
        int num_static_contacts = 0;

        for (int iaw : awake) {
            sheet.stone[iaw].set_goal_info();
            collidable[iaw] = (awake & ~BitSet32::from_index(iaw)) | asleep;
        }
        for (int ias0 : asleep) {
            for (int ias1 : BitSet32(asleep & ((1U << ias0) - 1))) {
                if (has_contact(sheet.stone[ias0], sheet.stone[ias1])) contact[num_static_contacts++] =
                    std::make_tuple(ias0, ias1, distance(sheet.stone[ias0], sheet.stone[ias1]), fpn_t(0));
            }
        }

        for (int t_ = 0; t_ < 1000; t_++) {
            // start collision procedure
            // ここから衝突処理
            int num_contacts = num_static_contacts;
            for (int iaw0 : awake) {
                const Stone& aw0 = sheet.stone[iaw0];
                for (int iaw1 : BitSet32(awake & collidable[iaw0] & ((1U << iaw0) - 1))) {
                    const Stone& aw1 = sheet.stone[iaw1];
                    if (has_contact(aw0, aw1)) contact[num_contacts++] =
                        std::make_tuple(iaw0, iaw1, distance(aw0, aw1), heading_speed(aw0, aw1));
                }
                for (int ias : BitSet32(asleep & collidable[iaw0])) {
                    const Stone& as = sheet.stone[ias];
                    if (has_contact(aw0, as)) contact[num_contacts++] =
                        std::make_tuple(iaw0, ias, distance(aw0, as), heading_speed(aw0, as, true));
                }
            }

            // if there is any stone contact involving a moving stone then collision process is done
            // 動いている石が関わる石の接触があれば衝突処理を行う
            if (num_contacts > num_static_contacts) {
                BitSet32 updated = resolve_contacts_(sheet, contact, num_contacts);

                awake = awake | updated;
                asleep = asleep & ~updated;

                // set collidable flag for every stone except stone that is just collided
                // 今衝突した石以外のすべての石に衝突可能フラグを立てる
                for (int iaw : updated) {
                    sheet.stone[iaw].set_goal_info();
                    collidable[iaw] = (awake & ~BitSet32::from_index(iaw)) | asleep;
                    for (int i = 0; i < num_contacts; i++) {
                        auto [i0, i1, _, __] = contact[i];
                        if      (i0 == iaw) collidable[iaw].reset(i1);
                        else if (i1 == iaw) collidable[iaw].reset(i0);
                    }
                }

                for (int iaw : BitSet32(awake & ~updated)) collidable[iaw] |= updated;

                // update static contacts
                // 静止コンタクトリストの整理
                for (int i = 0; i < num_static_contacts;) {
                    auto [i0, i1, _, __] = contact[i];
                    if (updated.test(i0) || updated.test(i1)) {
                        contact[i] = contact[--num_static_contacts];
                    } else i++;
                }
            }

            // calcurate lower bound time there can be no collision
            // 確実に衝突せず移動可能な時間を計算
            fpn_t min_contact_time = 100000;
            for (int iaw0 : awake) {
                const Stone& aw = sheet.stone[iaw0];
                fpn_t r_asleep = 100000;

                for (int ias : BitSet32(asleep & collidable[iaw0])) {
                    const Stone& as = sheet.stone[ias];
                    if (has_collision_chance_aw_as(aw, as)) {
                        r_asleep = std::min(r_asleep, min_contact_mr_aw_as(aw, as));
                    } else { // no collision possibility 非衝突が証明
                        collidable[iaw0].reset(ias);
                    }
                }
                if (r_asleep != 100000) {
                    r_asleep = std::min(r_asleep, aw.gr * (aw.gth - aw.th).c);
                    fpn_t ngr = std::sqrt(std::max(fpn_t(0), aw.gr * aw.gr + r_asleep * r_asleep - 2 * aw.gr * r_asleep * (aw.gth - aw.th).c));
                    fpn_t ngt = __r2t(ngr);
                    min_contact_time = std::min(min_contact_time, aw.gt - ngt);
                }

                for (int iaw1 : BitSet32(awake & collidable[iaw0] & ((1U << iaw0) - 1))) {
                    const Stone& aw1 = sheet.stone[iaw1];
                    fpn_t dt = min_contact_time_2aw(aw, aw1, min_contact_time);
                    if (dt == 100000) { // no collision possibility 非衝突が証明
                        collidable[iaw0].reset(iaw1);
                        collidable[iaw1].reset(iaw0);
                    } else {
                        min_contact_time = std::min(min_contact_time, dt);
                    }
                }
            }

            // if no stones are possible to collide, locate every moving stone to goal positions and finish
            // 全ての石が衝突せずに静止する事が証明されたならば、すべての石を置いて終わり
            if (min_contact_time == 100000) break;
            if (min_contact_time <= 0) break;

            // update each stone status by time
            // それぞれの石を移動可能な時間分動かす
            for (int iaw : awake) {
                Stone& aw = sheet.stone[iaw];
                if (aw.gt <= min_contact_time) {
                    // stop a stone whose stopping time is earlier than min_contact_time
                    // 最短衝突時刻より停止時刻が早い石は静止
                    aw.step_to_stop();
                    awake.reset(iaw);

                    if (!rink_only || is_in_play_area(aw.x, aw.y)) {
                        asleep.set(iaw);
                        for (int ias : asleep) {
                            const Stone& as = sheet.stone[ias];
                            if (has_contact(aw, as)) contact[num_static_contacts++] =
                                std::make_tuple(iaw, ias, distance(aw, as), fpn_t(0));
                        }
                    }
                } else {
                    aw.step_by_next_time(min_contact_time);
                }
            }
        }

        // locate stones and finish simulation
        // 停止位置に石を置いて終了
        for (int iaw : awake) {
            Stone& aw = sheet.stone[iaw];
            aw.step_to_stop();
            if (!rink_only || is_in_play_area(aw.x, aw.y)) asleep.set(iaw);
        }

        sheet.stone_bits = asleep;
        sheet.moving_stone_bits = 0;
        return false;
    }

    bool simulate(Sheet& sheet, int turn, const Shot& shot, bool freeguard, bool rink_only) {
        if (!sheet.stone_bits) { // free draw with no obstacles
            Position p = shot2dest(shot);
            if (!rink_only || is_in_play_area(p.x, p.y)) sheet.set_stone(turn, p);
            return false;
        } else {
            BitSet32 fg = 0;
            Sheet tmp;
            if (freeguard) {
                tmp = sheet;
                for (int i = 1 - turn % 2; i < 16; i += 2) {
                    if ((sheet.stone_bits & (1U << i)) && is_in_freeguard_zone(sheet.stone[i].x, sheet.stone[i].y)) fg.set(i);
                }
            }
            sheet.set_moving_stone(turn, X_THROW, Y_THROW, std::min(V_SIMULATION_MAX, xy2r(shot.vx, shot.vy)), xy2ang(shot.vx, shot.vy), spin2wthrow(shot.spin));
            simulate(sheet, rink_only);
            if (fg & ~sheet.stone_bits) {
                sheet = tmp; // unnecessary if original board is already copied before simulation
                return true;
            }
            return false;
        }
    }

    // table generation
    void init_table() {
        const double step = 0.0001;
        double ow = -5;
        double x = 0, y = 0, w = ow, vx = 0, vy = 5.01;
        int last_index = 0;
        int last_rest = 0;
        double last_time = 0;
        double last_vth = 0;
        std::vector<std::array<double, 6>> tmp(100001);
        for (int t = 0; t <= 1000000; t++) {
            if (!(vx == 0 && vy == 0))
            {
                last_vth = std::max(0.0, std::atan2(vx, vy));
                last_time = t * step;
                last_rest = t % 10;
            }
            if (t % 10 == 0) {
                int i = 100000 - t / 10;
                tmp[i][0] = x;
                tmp[i][1] = y;
                tmp[i][2] = std::sqrt(vx * vx + vy * vy);
                tmp[i][3] = last_vth;
                tmp[i][4] = w - ow;
                tmp[i][5] = last_time;
                if (vx == 0 && vy == 0) {
                    for (int j = 0; j < i; j++) {
                        for (int k = 0; k < 6; k++) tmp[j][k] = tmp[i][k];
                    }
                    last_index = i;
                    break;
                }
            }
            friction_step(vx, vy, w, &vx, &vy, &w, step / 2);
            x += vx * step;
            y += vy * step;
            friction_step(vx, vy, w, &vx, &vy, &w, step / 2);
        }
        //cerr << last_index << " " << last_vth << endl;
        for (int j = 0; j < last_index; j++) {
            for (int k = 0; k < 6; k++) tmp[j][k] = tmp[last_index][k];
        }

        // save to table
        // 到達位置を(0, 0)、到達時の向きをy軸負の向きとし回転と並行移動してテーブルに保存
        for (int i = 0; i <= 100000; i++) {
            double dx = x - tmp[i][0];
            double dy = y - tmp[i][1];
            double r = std::sqrt(dx * dx + dy * dy);
            double th = r == 0 ? 0 : (last_vth - std::atan2(dx, dy));
            double vth = last_vth - tmp[i][3];
            double t = last_time - tmp[i][5];
            assert(th >= -1e-4);
            assert(vth >= -1e-4);

            tmp[i][0] = r;
            tmp[i][1] = std::max(0.0, th);
            tmp[i][3] = std::max(0.0, vth);
            tmp[i][5] = t;
        }

        // table from time to stop
        // 時刻からの変換テーブル
        for (int i = 0; i <= __TABLE; i++) {
            int j = std::min(100000 - 1, last_index + i * (100000 / __TABLE));
            double c = i == 0 ? 0 : (1 - last_rest / 10.0);
            traj_table[i][0] = (1 - c) * tmp[j][0] + c * tmp[j + 1][0];
            traj_table[i][1] = (1 - c) * tmp[j][2] + c * tmp[j + 1][2];
            traj_table[i][2] = (1 - c) * tmp[j][4] + c * tmp[j + 1][4];
            double th  = (1 - c) * tmp[j][1] + c * tmp[j + 1][1];
            double vth = (1 - c) * tmp[j][3] + c * tmp[j + 1][3];
            traj_table[i][3] = std::sin(th);
            traj_table[i][4] = std::cos(th);
            traj_table[i][5] = std::sin(vth);
            traj_table[i][6] = std::cos(vth);
        }

        // table from velocity
        // 速度からの変換テーブル
        for (int i = 0; i <= __TABLE; i++) {
            double v = 5.0 * i / __TABLE;
            int j = i == 0 ? 0 : std::upper_bound(tmp.begin(), tmp.end() - 1, std::array<double, 6>({0, 0, v, 0, 0, 0}), [](auto x, auto y) { return x[2] < y[2]; }) - tmp.begin() - 1;
            double c = i == 0 ? 0.0 : (v - tmp[j][2]) / (tmp[j + 1][2] - tmp[j][2]);
            double r   = (1 - c) * tmp[j][0] + c * tmp[j + 1][0];
            double oth = (1 - c) * tmp[j][1] + c * tmp[j + 1][1];
            double vth = (1 - c) * tmp[j][3] + c * tmp[j + 1][3];
            double t   = (1 - c) * tmp[j][5] + c * tmp[j + 1][5];
            v2rtt_table[i][0] = r;
            v2rtt_table[i][1] = t;
            v2rtt_table[i][2] = std::sin(vth - oth);
            v2rtt_table[i][3] = std::cos(vth - oth);
            //cerr << v << " " << j << " " << c << " " << v2rtt_table[i][0] << " " << v2rtt_table[i][1] << " " << v2rtt_table[i][2] << endl;
        }

        // table from distance to stop
        // 距離からの変換テーブル
        for (int i = 0; i <= __TABLE; i++) {
            double r_ = 14.0 * i / __TABLE;
            double r = r_ * r_;
            int j = i == 0 ? 0 : std::upper_bound(tmp.begin(), tmp.end() - 1, std::array<double, 6>({r, 0, 0, 0, 0, 0}), [](auto x, auto y) { return x[0] < y[0]; }) - tmp.begin() - 1;
            double c = tmp[j][0] == tmp[j + 1][0] ? 0 : std::min(1.0, (r - tmp[j][0]) / (tmp[j + 1][0] - tmp[j][0]));
            double v = (1 - c) * tmp[j][2] + c * tmp[j + 1][2];
            double t = (1 - c) * tmp[j][5] + c * tmp[j + 1][5];
            r2vt_table[i][0] = v;
            r2vt_table[i][1] = t;
            //cerr << r << " " << j << " " << c << " " << r2vt_table[i][0] << " " << r2vt_table[i][1] << endl;
        }
    }
}
