#ifndef EXPANSION444_H_
#define EXPANSION444_H_

#include <octopus/vector4d.hpp>

#include <octopus/vector.h>

static const double delta[3][3] = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
static const size_t map2[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
static const size_t map3[3][3][3] = { { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } }, { { 1, 3, 4 }, { 3, 6, 7 }, { 4, 7, 8 } }, { { 2, 4, 5 }, { 4, 7, 8 }, { 5, 8,
        9 } } };

#define POLE0_CNT 1
#define POLE1_CNT (POLE0_CNT+1+2)
#define POLE2_CNT (POLE1_CNT+1+2+3)
#define POLE3_CNT (POLE2_CNT+1+2+3+4)

template<int EXPANSION_POLES, bool NO_DIPOLE = false>
class Expansion: public Vector<double, EXPANSION_POLES> {
    friend class boost::serialization::access;
    template<typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        for (boost::uint64_t i = 0; i < EXPANSION_POLES; i++) {
            ar & (*this)[i];
        }
    }

public:
    Expansion() {
    }

    double operator ()() const {
        return (*this)[0];
    }
    double& operator ()() {
        return (*this)[0];
    }

    double operator ()(int i) const {
        return (*this)[EXPANSION_POLES - 3 + i];
    }
    double& operator ()(int i) {
        return (*this)[EXPANSION_POLES - 3 + i];
    }

    double operator ()(int i, int j) const {
        return (*this)[1 + map2[i][j]];
    }
    double& operator ()(int i, int j) {
        return (*this)[1 + map2[i][j]];
    }

    double operator ()(int i, int j, int k) const {
        return (*this)[1 + 6 + map3[i][j][k]];
    }
    double& operator ()(int i, int j, int k) {
        return (*this)[1 + 6 + map3[i][j][k]];
    }

    Expansion& operator =(const Expansion& expansion) {
        for (int i = 0; i < EXPANSION_POLES; i++) {
            (*this)[i] = expansion[i];
        }
        return *this;
    }

    Expansion& operator =(double expansion) {
        for (int i = 0; i < EXPANSION_POLES; i++) {
            (*this)[i] = expansion;
        }
        return *this;
    }

    Expansion operator<<(const _3Vec& dX) const {
        Expansion you = *this;
        you <<= dX;
        return you;
    }

    Expansion& operator<<=(const _3Vec& dX) {
        Expansion& me = *this;
        for (int a = 0; a < 3; a++) {
            me() += me(a) * dX[a];
        }
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                me() += me(a, b) * dX[a] * dX[b] * 0.5;
            }
        }
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                me(a) += me(a, b) * dX[b];
            }
        }
        if (EXPANSION_POLES > POLE2_CNT) {
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    for (int c = 0; c < 3; c++) {
                        me() += me(a, b, c) * dX[a] * dX[b] * dX[c] * (1.0 / 6.0);
                    }
                }
            }
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    for (int c = 0; c < 3; c++) {
                        me(a) += me(a, b, c) * dX[b] * dX[c] * 0.5;
                    }
                }
            }
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    for (int c = a; c < 3; c++) {
                        me(a, c) += me(a, b, c) * dX[b];
                    }
                }
            }
        }
        return me;
    }

    Expansion operator>>(const _3Vec& dX) const {
        Expansion you = *this;
        you >>= dX;
        return you;
    }

    Expansion& operator>>=(const _3Vec& Y) {
        Expansion& me = *this;
        Expansion tmp;
        tmp = me;
        if (!NO_DIPOLE) {
            for (int p = 0; p < 3; p++) {
                me(p) += tmp() * Y[p];
            }
        }
        for (int p = 0; p < 3; p++) {
            for (int q = p; q < 3; q++) {
                if (!NO_DIPOLE) {
                    me(p, q) += tmp(q) * Y[p] + tmp(p) * Y[q];
                }
                me(p, q) += tmp() * Y[p] * Y[q];
            }
        }
        if (EXPANSION_POLES > POLE2_CNT) {
            for (int p = 0; p < 3; p++) {
                for (int q = p; q < 3; q++) {
                    for (int r = q; r < 3; r++) {
                        if (!NO_DIPOLE) {
                            me(p, q, r) += tmp(p) * Y[q] * Y[r] + tmp(q) * Y[r] * Y[p] + tmp(r) * Y[p] * Y[q];
                        }
                        me(p, q, r) += tmp() * Y[p] * Y[q] * Y[r] + tmp(p, q) * Y[r] + tmp(q, r) * Y[p] + tmp(r, p) * Y[q];
                    }
                }
            }
        }
        return me;
    }

    Vector<double, EXPANSION_POLES>& operator +=(const Vector<double, EXPANSION_POLES>& vec) {
        for (int i = 0; i < EXPANSION_POLES; i++) {
            (*this)[i] += vec[i];
        }
        return *this;
    }

    void compute_D(const _3Vec& Y) {
        Expansion& me = *this;
        const double r2inv = 1.0 / Y.dot(Y);
        const double d0 = -sqrt(r2inv);
        const double d1 = -d0 * r2inv;
        const double d2 = -3.0 * d1 * r2inv;
        const double d3 = -5.0 * d2 * r2inv;
        me() = d0;
        for (int a = 0; a < 3; a++) {
            me(a) = Y[a] * d1;
        }
        for (int a = 0; a < 3; a++) {
            for (int b = a; b < 3; b++) {
                me(a, b) = Y[a] * Y[b] * d2 + delta[a][b] * d1;
            }
        }
        for (int a = 0; a < 3; a++) {
            for (int b = a; b < 3; b++) {
                for (int c = b; c < 3; c++) {
                    me(a, b, c) = Y[a] * Y[b] * Y[c] * d3 + (delta[a][b] * Y[c] + delta[b][c] * Y[a] + delta[c][a] * Y[b]) * d2;
                }
            }
        }
    }

    void invert() {
        Expansion& me = *this;
        for (int a = 0; a < 3; a++) {
            me(a) = -me(a);
        }
        for (int a = 0; a < 3; a++) {
            for (int b = a; b < 3; b++) {
                for (int c = b; c < 3; c++) {
                    me(a, b, c) = -me(a, b, c);
                }
            }
        }
    }
    ~Expansion() {
    }
}
;

template<int EXPANSION_POLES>
inline Expansion<EXPANSION_POLES> operator*(const Expansion<EXPANSION_POLES>& M, const Expansion<EXPANSION_POLES>& D) {
    Expansion<EXPANSION_POLES> me;

    me() = M() * D();
    for (int a = 0; a < 3; a++) {
        me() -= M(a) * D(a);
    }
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            me() += M(a, b) * D(a, b) * 0.5;
        }
    }
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int c = 0; c < 3; c++) {
                me() -= M(a, b, c) * D(a, b, c) * (1.0 / 6.0);
            }
        }
    }

    for (int a = 0; a < 3; a++) {
        me(a) = M() * D(a);
    }
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            me(a) -= M(a) * D(a, b);
        }
    }
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int c = 0; c < 3; c++) {
                me(a) += M(c, b) * D(a, b, c) * 0.5;
            }
        }
    }

    for (int a = 0; a < 3; a++) {
        for (int b = a; b < 3; b++) {
            me(a, b) = M() * D(a, b);
        }
    }
    for (int a = 0; a < 3; a++) {
        for (int b = a; b < 3; b++) {
            for (int c = 0; c < 3; c++) {
                me(a, b) -= M(c) * D(a, b, c);
            }
        }
    }

    for (int a = 0; a < 3; a++) {
        for (int b = a; b < 3; b++) {
            for (int c = b; c < 3; c++) {
                me(a, b, c) = M() * D(a, b, c);
            }
        }
    }

    return me;
}

template<int EXPANSION_POLES>
inline Expansion<EXPANSION_POLES> operator*(const Expansion<EXPANSION_POLES - 3, true>& M, const Expansion<EXPANSION_POLES>& D) {
    Expansion<EXPANSION_POLES> me;

    me() = M() * D();
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            me() += M(a, b) * D(a, b) * 0.5;
        }
    }
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int c = 0; c < 3; c++) {
                me() -= M(a, b, c) * D(a, b, c) * (1.0 / 6.0);
            }
        }
    }

    for (int a = 0; a < 3; a++) {
        me(a) = M() * D(a);
    }
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            for (int c = 0; c < 3; c++) {
                me(a) += M(c, b) * D(a, b, c) * 0.5;
            }
        }
    }

    for (int a = 0; a < 3; a++) {
        for (int b = a; b < 3; b++) {
            me(a, b) = M() * D(a, b);
        }
    }

    for (int a = 0; a < 3; a++) {
        for (int b = a; b < 3; b++) {
            for (int c = b; c < 3; c++) {
                me(a, b, c) = M() * D(a, b, c);
            }
        }
    }

    return me;
}

typedef Expansion<POLE3_CNT> expansion_t;
typedef Expansion<POLE3_CNT - 3, true> expansion_nodip_t;
typedef Expansion<POLE2_CNT> expansion_m1_t;

struct multipole_t {
    expansion_nodip_t M;
    _3Vec X;
    bool is_leaf;

    friend class boost::serialization::access;

    template<typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & X[0] & X[1] & X[2] & M & is_leaf;
    }

};

struct taylor_t {
    expansion_t phi;
    expansion_m1_t g_lz;
    double x, y, z;
};

struct multipole_dot_t {
    expansion_t M_dot;
};

struct taylor_dot_t {
    expansion_t phi_dot;
};

typedef octopus::vector4d<multipole_t, 1> multipole_array_t;

#endif /* EXPANSION_H_ */

