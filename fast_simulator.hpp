
// Ayumu Curling Simulator
// Katsuki Ohto

#ifndef AYUMU_CURLING_SIMULATOR_HPP_
#define AYUMU_CURLING_SIMULATOR_HPP_

#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif


namespace FastSimulator {
    using fpn_t = double;

    static const char *C = "BW";

    // rink parameters
    const fpn_t X_THROW = 0;
    const fpn_t Y_THROW = 0;

    const fpn_t R_STONE_RADIUS = 0.145;
    const fpn_t R_HOUSE_RADIUS = 1.829;

    const fpn_t R_PA_WIDTH = 4.75;
    const fpn_t R_PA_LENGTH = 8.23;

    const fpn_t X_TEE = X_THROW;
    const fpn_t Y_TEE = Y_THROW + 38.405;

    const fpn_t X_PA_MAX = X_TEE + R_PA_WIDTH / 2;
    const fpn_t Y_PA_MAX = Y_TEE + R_HOUSE_RADIUS;
    const fpn_t X_PA_MIN = X_PA_MAX - R_PA_WIDTH;
    const fpn_t Y_PA_MIN = Y_PA_MAX - R_PA_LENGTH;

    const fpn_t X_RINK_MIN = X_PA_MIN;
    const fpn_t Y_RINK_MIN = Y_THROW - R_STONE_RADIUS * 2;
    const fpn_t X_RINK_MAX = X_PA_MAX;
    const fpn_t Y_RINK_MAX = Y_PA_MAX + R_STONE_RADIUS * 2;

    // stone parameters
    const fpn_t FRICTION_STONES = 0.2;

    const fpn_t RESTITUTION_STONES = 1; // perfectly elastic collision
    const fpn_t RESTITUTION_THRESHOLD = 0;

    const fpn_t MASS_STONE = 19.96;

    const fpn_t DENSITY_STONE = MASS_STONE / (M_PI * R_STONE_RADIUS * R_STONE_RADIUS);

    const fpn_t W_THROW = M_PI / 2;

    inline bool is_left_spin(int spin) { return spin != 0; }

    inline fpn_t spin2wthrow(int spin) { return is_left_spin(spin) ? W_THROW : -W_THROW; }

    template <typename float_t>
    int w2spin(float_t w) { return std::signbit(w) ? 0 : 1; } // w <= 0 -> right

    const fpn_t R_COLLISION = 0.00004; // distance threshold to determine collision

    const fpn_t V_SIMULATION_MAX = 4.9;

    // basic functions
    template <typename float_t>
    float_t r2r2(float_t r) {
        return r * r;
    }

    template <typename float_t>
    float_t xy2r2(float_t x, float_t y) {
        return x * x + y * y;
    }

    template <typename float_t>
    float_t xy2r(float_t x, float_t y) {
        return std::sqrt(x * x + y * y);
    }

    template <typename float_t>
    float_t xy2th(float_t x, float_t y) {
        return std::atan2(x, y);
    }

    template <typename T>
    T sign(T v) {
        return v > 0 ? T(1) : (v < 0 ? T(-1) : T(0));
    }

    template <typename T>
    void rotate(T ox, T oy, T theta, T *x, T *y) {
        *x = ox * std::cos(theta) - oy * std::sin(theta);
        *y = oy * std::cos(theta) + ox * std::sin(theta);
    }

    // basic structures
    template <typename T>
    struct MiniBitSet {
        T data;

        constexpr MiniBitSet(): data() {}
        constexpr MiniBitSet(T d): data(d) {}
        constexpr MiniBitSet(const MiniBitSet<T>& bs): data(bs.data) {}

        static constexpr MiniBitSet<T> from_index(int index) { return MiniBitSet<T>(T(1) << index); }

        constexpr bool empty() const { return !data; }
        constexpr operator T() const { return data; }

        constexpr bool test(int index) const { return (data >> index) & T(1); }
        int count() const { return std::popcount(data); }

        constexpr bool holds(const MiniBitSet<T>& bs) const { return !(~data & bs.data); }
        constexpr bool exclusive(const MiniBitSet<T>& bs) const { return !(data & bs.data); }

        MiniBitSet<T>& operator |=(T d) { data |= d; return *this; }
        MiniBitSet<T>& operator &=(T d) { data &= d; return *this; }
        MiniBitSet<T>& operator ^=(T d) { data ^= d; return *this; }
        MiniBitSet<T>& operator +=(T d) { data += d; return *this; }
        MiniBitSet<T>& operator -=(T d) { data -= d; return *this; }
        MiniBitSet<T>& operator <<=(int i) { data <<= i; return *this; }
        MiniBitSet<T>& operator >>=(int i) { data >>= i; return *this; }
        MiniBitSet<T>& operator <<=(unsigned int i) { data <<= i; return *this; }
        MiniBitSet<T>& operator >>=(unsigned int i) { data >>= i; return *this; }

        MiniBitSet<T>& reset() { data = 0; return *this; }
        MiniBitSet<T>& set() { data = T(-1); return *this; }
        MiniBitSet<T>& reset(int index) { data &= ~(T(1) << index); return *this; }
        MiniBitSet<T>& set(int index) { data |= T(1) << index; return *this; }

        int lowest() const { assert(!empty()); return std::countr_zero(data); }
        int highest() const { assert(!empty()); return sizeof(T) * 8 - 1 - std::countl_zero(data); }
        int pop_lowest() { int index = lowest(); data &= data - T(1); return index; }
        int pop_highest() { int index = highest(); reset(index); return index; }

        class const_iterator {
            friend MiniBitSet<T>;
        public:
            using difference_type   = std::ptrdiff_t;
            using value_type        = int;
            using pointer           = int*;
            using reference         = int&;
            using iterator_category = std::input_iterator_tag;
            int operator *() const { return std::countr_zero(data_); }
            bool operator !=(const const_iterator& itr) const { return pclass_ != itr.pclass_ || data_ != itr.data_; }
            const_iterator& operator ++() { data_ &= data_ - T(1); return *this; }
        protected:
            explicit const_iterator(const MiniBitSet<T> *pclass): pclass_(pclass), data_(pclass->data) {}
            explicit const_iterator(const MiniBitSet<T> *pclass, T d): pclass_(pclass), data_(d) {}
            const MiniBitSet<T> *const pclass_;
            T data_;
        };

        const_iterator begin() const { return const_iterator(this); }
        const_iterator end() const { return const_iterator(this, 0); }
    };

    template <typename T>
    std::ostream& operator <<(std::ostream& ost, const MiniBitSet<T>& bs) {
        ost << "{";
        int count = 0;
        for (int index : bs) { ost << (count ? " " : "") << index; count += 1; }
        ost << "}";
        return ost;
    }

    using BitSet32 = MiniBitSet<std::uint32_t>;

    // basic functions for curling
    template <typename float_t>
    bool is_on_rink(float_t x, float_t y) {
        return X_RINK_MIN + R_STONE_RADIUS < x && x < X_RINK_MAX - R_STONE_RADIUS
            && Y_RINK_MIN + R_STONE_RADIUS < y && y < Y_RINK_MAX - R_STONE_RADIUS;
    }

    template <typename float_t>
    bool is_in_play_area(float_t x, float_t y) {
        return X_PA_MIN + R_STONE_RADIUS < x && x < X_PA_MAX - R_STONE_RADIUS
            && Y_PA_MIN + R_STONE_RADIUS < y && y < Y_PA_MAX + R_STONE_RADIUS;
    }

    template <typename float_t>
    bool is_in_house(float_t x, float_t y) {
        return xy2r2(x - X_TEE, y - Y_TEE) < r2r2(R_HOUSE_RADIUS + R_STONE_RADIUS);
    }

    template <typename float_t>
    bool is_in_freeguard_zone(float_t x, float_t y) {
        if (!(X_PA_MIN + R_STONE_RADIUS < x && x < X_PA_MAX - R_STONE_RADIUS &&
              Y_PA_MIN + R_STONE_RADIUS < y && y < Y_TEE - R_STONE_RADIUS)) return false;
        if (is_in_house(x, y)) return false;
        return true;
    }

    // structures for curling
    struct Position {
        fpn_t x, y;
        Position() {}
        Position(fpn_t x_, fpn_t y_): x(x_), y(y_) {}
    };

    struct PositionTime : public Position {
        fpn_t t;
        PositionTime() {}
        PositionTime(fpn_t x, fpn_t y, fpn_t t_): Position(x, y), t(t_) {}
    };

    inline std::ostream& operator <<(std::ostream& ost, const Position& pos) {
        ost << "(" << pos.x << ", " << pos.y << ")";
        return ost;
    }

    inline fpn_t distance2(const Position& p0, const Position& p1) {
        return xy2r2(p1.x - p0.x, p1.y - p0.y);
    }

    inline fpn_t distance(const Position& p0, const Position& p1) {
        return xy2r(p1.x - p0.x, p1.y - p0.y);
    }

    inline Position random_position_in_house(fpn_t r0, fpn_t r1) {
        fpn_t r = (R_HOUSE_RADIUS + R_STONE_RADIUS) * r0;
        return Position(
            X_TEE + r * std::sin(2 * M_PI * r1),
            Y_TEE + r * std::cos(2 * M_PI * r1)
        );
    }

    inline Position random_position_in_play_area(fpn_t r0, fpn_t r1) {
        return Position(
            X_PA_MIN + R_STONE_RADIUS + (X_PA_MAX - X_PA_MIN - 2 * R_STONE_RADIUS) * r0,
            Y_PA_MIN + R_STONE_RADIUS + (Y_PA_MAX - Y_PA_MIN) * r1
        );
    }

    inline bool has_contact(const Position& p0, const Position& p1) {
        return distance2(p0, p1) < r2r2(R_STONE_RADIUS * 2 + R_COLLISION);
    }

    inline bool has_bad_contact(const Position& p0, const Position& p1) {
        return distance2(p0, p1) < r2r2(R_STONE_RADIUS * 2);
    }

    inline fpn_t stone_gap(const Position& p0, const Position& p1) {
        return distance(p0, p1) - 2 * R_STONE_RADIUS;
    }

    struct Shot {
        fpn_t vx, vy;
        int spin;
        Shot() {}
        Shot(fpn_t vx_, fpn_t vy_, int spin_): vx(vx_), vy(vy_), spin(spin_) {}
        static Shot from_vth(fpn_t v, fpn_t vth, int spin_) {
            return Shot(v * std::sin(vth), v * std::cos(vth), spin_);
        }
    };

    inline std::ostream& operator <<(std::ostream& ost, const Shot& shot) {
        ost << "(" << shot.vx << ", " << shot.vy << ", " << shot.spin << ")";
        return ost;
    }

    struct Angle {
        fpn_t s, c;
        Angle() {}
        Angle(fpn_t s_, fpn_t c_): s(s_), c(c_) {}
        Angle(fpn_t th): s(std::sin(th)), c(std::cos(th)) {}
        static Angle from_cos(fpn_t c_) { return Angle(std::sqrt(1 - c_ * c_), c_); }
        fpn_t radian() const { return std::atan2(s, c); }
    };

    inline std::ostream& operator <<(std::ostream& ost, const Angle& a) {
        ost << a.radian();
        return ost;
    }

    inline Angle operator -(const Angle& a) { return Angle(-a.s, a.c); }
    inline Angle operator +(const Angle& a, const Angle& b) { return Angle(a.s * b.c + a.c * b.s, a.c * b.c - a.s * b.s); }
    inline Angle operator -(const Angle& a, const Angle& b) { return Angle(a.s * b.c - a.c * b.s, a.c * b.c + a.s * b.s); }
    inline bool operator <(const Angle& a, const Angle& b) { return (b - a).s > 0; }
    inline bool operator >(const Angle& a, const Angle& b) { return (b - a).s < 0; }

    inline Angle xy2ang(fpn_t x, fpn_t y) { fpn_t r = xy2r(x, y); return Angle(x / r, y / r); }

    inline void rotate(fpn_t ox, fpn_t oy, const Angle& theta, fpn_t *x, fpn_t *y) {
        *x = ox * theta.c - oy * theta.s;
        *y = oy * theta.c + ox * theta.s;
    }

    struct Stone : public Position {
        fpn_t v, w;
        Angle th;

        fpn_t gx, gy;
        fpn_t gr, gt;
        fpn_t wp;
        Angle gth, oth;

        fpn_t vx() const { return v * th.s; }
        fpn_t vy() const { return v * th.c; }

        void set(const Position& pos) {
            x = pos.x;
            y = pos.y;
            v = w = 0;
            th = Angle(0, 0);
        }

        void set(fpn_t x_, fpn_t y_, fpn_t v_, Angle th_, fpn_t w_) {
            x = x_;
            y = y_;
            v = v_;
            th = th_;
            w = w_;
        }

        std::string goal_info() const {
            std::ostringstream oss;
            oss << "(gx " << gx << " gy " << gy << " gr " << gr << " gth " << gth << " gt " << gt;
            oss << " oth " << oth << " wp " << wp << ")";
            return oss.str();
        }

        // implement after tables
        void set_goal_info();
        void step_to_stop();
        void step_by_next_time(fpn_t nt);
    };

    inline std::ostream& operator <<(std::ostream& ost, const Stone& s) {
        ost << "(" << s.x << ", " << s.y << "; " << s.v << ", " << s.th << ", " << s.w << ")";
        return ost;
    }

    inline bool is_heading(const Stone& s0, const Stone& s1, bool stopped = false) {
        fpn_t dx = s1.x - s0.x;
        fpn_t dy = s1.y - s0.y;
        fpn_t dvx = (stopped ? 0 : s1.vx()) - s0.vx();
        fpn_t dvy = (stopped ? 0 : s1.vy()) - s0.vy();
        return dx * dvx + dy * dvy < 0;
    }

    inline fpn_t heading_speed(const Stone& s0, const Stone& s1, bool stopped = false) {
        Angle coltheta = xy2ang(s1.x - s0.x, s1.y - s0.y);
        fpn_t v0n = s0.v * (s0.th - coltheta).c;
        fpn_t v1n = stopped ? 0 : s1.v * (s1.th - coltheta).c;
        return v0n - v1n;
    }

    struct Sheet {
        std::uint32_t stone_bits = 0;
        std::uint32_t moving_stone_bits = 0;

        // 0-th is the first stone of each end, and 0,2,4,7,8,10,12,14 is first team's stones
        // (this is important for the free guard rule)
        // 0番目をエンド最初の石とし、0,2,4,8,10,12,14が先手番の石と想定（フリーガードルールで重要）
        Stone stone[16];

        void clear() { stone_bits = moving_stone_bits = 0; }

        template <class position_t>
        void set_stone(int index, const position_t& pos) {
            stone[index].set(pos);
            stone_bits |= 1U << index;
        }
        void set_moving_stone(int index, fpn_t x, fpn_t y, fpn_t v, Angle th, fpn_t w) {
            stone[index].set(x, y, v, th, w);
            moving_stone_bits |= 1U << index;
            stone_bits |= 1U << index;
        }
    };

    inline std::ostream& operator <<(std::ostream& ost, const Sheet& sheet) {
        for (int c = 0; c < 2; c++) {
            ost << C[c];
            for (int i = 0; i < 8; i++) {
                int index = i * 2 + c;
                if (sheet.stone_bits & (1U << index)) {
                    ost << index << Position(sheet.stone[index]);
                }
            }
        }
        return ost;
    }

    inline bool placable(const Sheet& sheet, const Position& pos) {
        for (int index : BitSet32(sheet.stone_bits)) {
            if (has_bad_contact(sheet.stone[index], pos)) return false;
        }
        return true;
    }

    // simulation base
    template <typename float_t>
    void friction_step(float_t vx, float_t vy, float_t w, float_t *nvx, float_t *nvy, float_t *nw, float_t time_step) {
        // update velocity
        float_t v = xy2r(vx, vy);
        float_t dv = (float_t(0.00200985) / (v + float_t(0.06385782)) + float_t(0.00626286)) * float_t(9.80665) * time_step;
        float_t new_v = std::max(float_t(0), v - dv);

        if (new_v == 0) {
            *nvx = *nvy = *nw = 0;
        } else {
            float_t yaw = sign(w) * float_t(0.00820) * std::pow(v, float_t(-0.8)) * time_step;
            float_t v_long = new_v * std::cos(yaw);
            float_t v_trans = new_v * std::sin(yaw);
            float_t e_long[2] = {vx / v, vy / v};
            float_t e_trans[2] = {-e_long[1], e_long[0]};
            *nvx = v_long * e_long[0] + v_trans * e_trans[0];
            *nvy = v_long * e_long[1] + v_trans * e_trans[1];

            float_t dw = float_t(0.025) / std::max(v, float_t(0.001)) * time_step;
            *nw = std::abs(w) <= dw ? 0 : w - dw * sign(w);
        }
    }

    extern PositionTime simulate_solo(Shot shot, fpn_t time_step);

    // collision
    extern void collision(Stone *s0, Stone *s1, bool stopped = false);
    extern BitSet32 resolve_contacts(Sheet& sheet);

    // collision possibility
    extern bool has_collision_chance_aw_as(const Stone& aw, const Stone& as);
    extern fpn_t min_contact_mr_aw_as(const Stone& aw, const Stone& as);
    extern fpn_t min_contact_time_2aw(const Stone& oaw0, const Stone& oaw1, fpn_t basetime);

    // trajectory data
    constexpr int __TABLE = 25000; // divisors of 100,000 only 100000の約数である必要あり
    extern float traj_table[__TABLE + 1][7]; // 0 <= t <= 100
    extern float v2rtt_table[__TABLE + 1][4]; // 0 <= v <= 5
    extern float r2vt_table[__TABLE + 1][2]; // 0 <= r <= 196

    // using trajectory
    inline fpn_t __load_traj_table(int i, fpn_t c, int k) {
        return (1 - c) * traj_table[i][k] + c * traj_table[i + 1][k];
    }
    inline fpn_t __t2_traj_table(fpn_t t, int k) {
        int i = t / fpn_t(100.0 / __TABLE);
        fpn_t c = t / fpn_t(100.0 / __TABLE) - i;
        return __load_traj_table(i, c, k);
    }

    inline fpn_t __t2r(fpn_t t) { return __t2_traj_table(t, 0); }
    inline fpn_t __t2v(fpn_t t) { return __t2_traj_table(t, 1); }
    inline fpn_t __t2wp(fpn_t t) { return __t2_traj_table(t, 2); }
    inline Angle __t2gth(fpn_t t) { return Angle(__t2_traj_table(t, 3), __t2_traj_table(t, 4)); }
    inline Angle __t2vth(fpn_t t) { return Angle(__t2_traj_table(t, 5), __t2_traj_table(t, 6)); }

    inline fpn_t __load_v2rtt_table(int i, fpn_t c, int k) {
        return (1 - c) * v2rtt_table[i][k] + c * v2rtt_table[i + 1][k];
    }
    inline fpn_t __v2_v2rtt_table(fpn_t v, int k) {
        int i = v / fpn_t(5.0 / __TABLE);
        fpn_t c = v / fpn_t(5.0 / __TABLE) - i;
        return __load_v2rtt_table(i, c, k);
    }

    inline fpn_t __v2r(fpn_t v) { return __v2_v2rtt_table(v, 0); }
    inline fpn_t __v2t(fpn_t v) { return __v2_v2rtt_table(v, 1); }
    inline Angle __v2th(fpn_t v) { return Angle(__v2_v2rtt_table(v, 2), __v2_v2rtt_table(v, 3)); }

    inline fpn_t __load_r2vt_table(int i, fpn_t c, int k) {
        return (1 - c) * r2vt_table[i][k] + c * r2vt_table[i + 1][k];
    }
    inline fpn_t __r2_r2vt_table(fpn_t r, int k) {
        r = std::sqrt(r);
        int i = r / fpn_t(14.0 / __TABLE);
        fpn_t c = r / fpn_t(14.0 / __TABLE) - i;
        return __load_r2vt_table(i, c, k);
    }

    inline fpn_t __r2v(fpn_t r) { return __r2_r2vt_table(r, 0); }
    inline fpn_t __r2t(fpn_t r) { return __r2_r2vt_table(r, 1); }

    inline Angle __vr2th(fpn_t v, fpn_t r) {
        fpn_t gt = __v2t(v);
        fpn_t gr = __t2r(gt);
        Angle gth = __t2gth(gt);
        Angle gvth = __t2vth(gt);
        fpn_t dr = gr - r;
        while (1) {
            fpn_t dt = __r2t(dr);
            Angle dth = __t2gth(dt);
            fpn_t nr2 = gr * gr + dr * dr - 2 * gr * dr * (gth - dth).c;
            fpn_t nr = std::sqrt(std::max(fpn_t(0), nr2));
            dr += nr - r;
            if (std::abs(nr - r) < 1e-4) break;
        }
        Angle th0 = gvth - gth;
        Angle th1 = Angle::from_cos((gr * gr + r * r - dr * dr) / (2 * gr * r));
        return th0 - th1;
    }

    // basic simulation functions
    inline Position shot2dest(const Shot& shot) {
        // calculate stopping position of a specified shot
        // ショットの停止位置の計算
        fpn_t v = xy2r(shot.vx, shot.vy);
        Angle vtheta = xy2ang(shot.vx, shot.vy);
        fpn_t r = __v2r(v);
        Angle dtheta = __v2th(v);
        if (is_left_spin(shot.spin)) dtheta = -dtheta;
        fpn_t dx, dy;
        rotate(r, fpn_t(0), vtheta + dtheta, &dy, &dx);
        return Position(X_THROW + dx, Y_THROW + dy);
    }

    inline Shot dest2shot(const Position& pos, int spin) {
        // calculate velocity (Vx, Vy) that stops at a specified position
        // the basic function for draw shots
        // 指定した位置で止まるショットの速度ベクトルの計算
        // ドロー系ショットの基本関数
        fpn_t dx = pos.x - X_THROW;
        fpn_t dy = pos.y - Y_THROW;
        fpn_t r = xy2r(dx, dy);
        Angle dtheta = xy2ang(dx, dy);
        fpn_t v = __r2v(r);
        Angle vtheta = __v2th(v);
        if (is_left_spin(spin)) vtheta = -vtheta;
        fpn_t vx, vy;
        rotate(v, fpn_t(0), dtheta - vtheta, &vy, &vx);
        return Shot(vx, vy, spin);
    }

    inline Shot passpoint2shot(const Position& pos, fpn_t v, int spin) {
        // calculate velocity (Vx, Vy) that passes through a specified position
        // with a specified initial speed
        // called in hit-type shots generation, mainly for take-out shots
        // 指定した初速で指定した位置を通るショットの速度ベクトルを計算
        // ヒット系、主にテイクアウト系の着手の生成で利用
        fpn_t dx = pos.x - X_THROW;
        fpn_t dy = pos.y - Y_THROW;
        fpn_t r = xy2r(dx, dy);
        Angle dtheta = xy2ang(dx, dy);
        Angle theta = __vr2th(v, r);
        if (is_left_spin(spin)) theta = -theta;
        fpn_t vx, vy;
        rotate(v, fpn_t(0), dtheta - theta, &vy, &vx);
        return Shot(vx, vy, spin);
    }

    inline Shot passpointgo2shot(const Position& pos, fpn_t dr, int spin) {
        // calculate velocity (Vx, Vy) that passes through a specified position
        // and then advances a specified distance (from the point of throw).
        // mainly called in pushing shots generation
        // 指定された位置を通り、さらに投げた点から見て指定距離分進む速度ベクトルを計算
        // 主に石を押し込むヒット系の着手の生成で呼ばれる
        fpn_t dx = pos.x - X_THROW;
        fpn_t dy = pos.y - Y_THROW;
        fpn_t r = xy2r(dx, dy);
        fpn_t v = __r2v(r + dr);
        Angle dtheta = xy2ang(dx, dy);
        Angle theta = __vr2th(v, r);
        if (is_left_spin(spin)) theta = -theta;
        fpn_t vx, vy;
        rotate(v, fpn_t(0), dtheta - theta, &vy, &vx);
        return Shot(vx, vy, spin);
    }

    // fast simulation
    extern bool simulate(Sheet& sheet, bool rink_only = true);
    extern bool simulate(Sheet& sheet, int turn, const Shot& shot, bool freeguard = false, bool rink_only = true);

    // table generation
    void init_table();

    struct Initializer {
        Initializer() {
            init_table();
        }
    };

    extern Initializer __initializer;
}

#endif // AYUMU_CURLING_SIMULATOR_HPP_
