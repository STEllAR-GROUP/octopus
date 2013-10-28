#ifndef VECTOR_H_
#define VECTOR_H_

#include "assert.h"

template<class T, int N>
class Vector {
private:
	T a[N];
public:
	Vector() {
		return;
	}
	~Vector() {
		return;
	}
	int size() const {
		return N;
	}
	bool operator==(const Vector<T, N>& v) const {
		bool rc = true;

		for (int i = 0; i < N; i++) {
			if (v.a[i] != a[i]) {
				rc = false;
				break;
			}
		}
		return rc;
	}
	bool operator!=(const Vector<T, N>& v) const {
		return !(v == *this);
	}
	Vector(const Vector<T, N>& v) {
		*this = v;
	}
	Vector(const T& v) {
		*this = v;
	}
	Vector<T, N>& operator=(const Vector<T, N>& v) {
		for (int i = 0; i < N; i++) {
			a[i] = v.a[i];
		}
		return *this;
	}
	T dot(const Vector<T, N>& v) const {
		T d = 0;
		for (int i = 0; i < N; i++) {
			d += a[i] * v.a[i];
		}
		return d;
	}
	T mag() const {
		return sqrt(this->dot(*this));
	}
	Vector<T, N>& operator=(const T& v) {
		for (int i = 0; i < N; i++) {
			a[i] = v;
		}
		return *this;
	}
	Vector<T, N>& operator/=(const T& t) {

		for (int i = 0; i < N; i++) {
			a[i] /= t;
		}
		return *this;
	}
	Vector<T, N>& operator*=(const T& t) {

		for (int i = 0; i < N; i++) {
			a[i] *= t;
		}
		return *this;
	}
	Vector<T, N>& operator-=(const T& t) {

		for (int i = 0; i < N; i++) {
			a[i] -= t;
		}
		return *this;
	}
	Vector<T, N>& operator+=(const T& t) {

		for (int i = 0; i < N; i++) {
			a[i] += t;
		}
		return *this;
	}
	Vector<T, N>& operator+=(const Vector<T, N>& v) {

		for (int i = 0; i < N; i++) {
			a[i] += v.a[i];
		}
		return *this;
	}
	Vector<T, N>& operator-=(const Vector<T, N>& v) {

		for (int i = 0; i < N; i++) {
			a[i] -= v.a[i];
		}
		return *this;
	}
	Vector<T, N>& operator*=(const Vector<T, N>& v) {

		for (int i = 0; i < N; i++) {
			a[i] *= v.a[i];
		}
		return *this;
	}
	Vector<T, N>& operator/=(const Vector<T, N>& v) {

		for (int i = 0; i < N; i++) {
			a[i] /= v.a[i];
		}
		return *this;
	}
	const Vector<T, N> operator+(const Vector<T, N>& v) const {
		return (Vector<T, N>(*this) += v);
	}
	const Vector<T, N> operator-(const Vector<T, N>& v) const {
		return (Vector<T, N>(*this) -= v);
	}
	const Vector<T, N> operator*(const Vector<T, N>& v) const {
		return (Vector<T, N>(*this) *= v);
	}
	const Vector<T, N> operator/(const Vector<T, N>& v) const {
		return (Vector<T, N>(*this) /= v);
	}
	const Vector<T, N> operator+(const T& v) const {
		return (Vector<T, N>(*this) += v);
	}
	const Vector<T, N> operator-(const T& v) const {
		return (Vector<T, N>(*this) -= v);
	}
	const Vector<T, N> operator*(const T& v) const {
		return (Vector<T, N>(*this) *= v);
	}
	const Vector<T, N> operator/(const T& v) const {
		return (Vector<T, N>(*this) /= v);
	}
	Vector<T, N> operator-() const {
		return (Vector<T, N>(T(0)) -= *this);
	}
	Vector<T, N> operator+() const {
		return Vector<T, N>(*this);
	}
	T& operator[](int i) {
		assert( i < N);
		assert( i >= 0);
		return a[i];
	}
	const T& operator[](int i) const {
		assert( i < N);
		assert( i >= 0);
		return a[i];
	}
	Vector<T, N>& operator <<=(int i) {

		for (int j = 0; j < N; j++) {
			a[j] <<= i;
		}
		return *this;
	}
	Vector<T, N>& operator >>=(int i) {

		for (int j = 0; j < N; j++) {
			a[j] >>= i;
		}
		return *this;
	}
};

template<class T, int N>
inline Vector<T, N> minmod(const Vector<T, N>& v1, const Vector<T, N>& v2) {
	Vector<T, N> mm;

	for (int i = 0; i < N; i++) {
		mm[i] = minmod(v1[i], v2[i]);
	}
	return mm;
}

template<class T, int N>
inline Vector<T, N> minmod(const Vector<T, N>& v1, const Vector<T, N>& v2, const Vector<T, N>& v3) {
	Vector<T, N> mm;

	for (int i = 0; i < N; i++) {
		mm[i] = minmod(v1[i], v2[i], v3[i]);
	}
	return mm;
}

template<class T, int N>
inline Vector<T, N> minmod_theta(const Vector<T, N>& v1, const Vector<T, N>& v2, T theta) {
	Vector<T, N> mm;

	for (int i = 0; i < N; i++) {
		mm[i] = minmod(theta * v1[i], theta * v2[i], 0.5 * (v1[i] + v2[i]));
	}
	return mm;
}

template<class T, int N>
Vector<int, N> nint(const Vector<T, N>& a) {
	Vector<int, N> v;

	for (int i = 0; i < N; i++) {
		v[i] = nint(a[i]);
	}
	return v;
}

typedef Vector<double, 3> _3Vec;

#endif /* VECTOR_H_ */
