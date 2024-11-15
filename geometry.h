#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <cmath>
#include <iostream>
#include <cassert>

template <class t> struct Vec2 {
    t x, y;
    Vec2<t>() : x(t()), y(t()) {}
    Vec2<t>(t _x, t _y) : x(_x), y(_y) {}
    Vec2<t>(const Vec2<t> &v) : x(t()), y(t()) { *this = v; }
    Vec2<t> & operator =(const Vec2<t> &v) {
        if (this != &v) {
            x = v.x;
            y = v.y;
        }
        return *this;
    }
    Vec2<t> operator +(const Vec2<t> &V) const { return Vec2<t>(x+V.x, y+V.y); }
    Vec2<t> operator -(const Vec2<t> &V) const { return Vec2<t>(x-V.x, y-V.y); }
    Vec2<t> operator *(float f)          const { return Vec2<t>(x*f, y*f); }
    t& operator[](const int i) {
      if (i == 0)
        return x;
      else
        return y;
    }
    template <class > friend std::ostream& operator<<(std::ostream& s, Vec2<t>& v);
};

template <class t> struct Vec3 {
    t x, y, z;
    Vec3<t>() : x(t()), y(t()), z(t()) { }
    Vec3<t>(t _x, t _y, t _z) : x(_x), y(_y), z(_z) {}
    template <class u> Vec3<t>(const Vec3<u> &v);
    Vec3<t>(const Vec3<t> &v) : x(t()), y(t()), z(t()) { *this = v; }
    Vec3<t> & operator =(const Vec3<t> &v) {
        if (this != &v) {
            x = v.x;
            y = v.y;
            z = v.z;
        }
        return *this;
    }
    Vec3<t> operator ^(const Vec3<t> &v) const { return Vec3<t>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); }
    Vec3<t> operator +(const Vec3<t> &v) const { return Vec3<t>(x+v.x, y+v.y, z+v.z); }
    Vec3<t> operator -(const Vec3<t> &v) const { return Vec3<t>(x-v.x, y-v.y, z-v.z); }
    Vec3<t> operator *(float f)          const { return Vec3<t>(x*f, y*f, z*f); }
    t       operator *(const Vec3<t> &v) const { return x*v.x + y*v.y + z*v.z; }
    float norm () const { return std::sqrt(x*x+y*y+z*z); }
    Vec3<t> & normalize(t l=1) { *this = (*this)*(l/norm()); return *this; }
    t& operator[](const int i) { if (i<=0) return x; else if (i==1) return y; else return z; }
    template <class > friend std::ostream& operator<<(std::ostream& s, Vec3<t>& v);
};

template <class t>
struct Vec4 {
  t x, y, z, w;
  Vec4<t>() : x(t()), y(t()), z(t()), w(t()) {}
  Vec4<t>(t _x, t _y, t _z, t _w) : x(_x), y(_y), z(_z), w(_w) {}
  Vec4<t>(const Vec4<t>& v) : x(v.x), y(v.y), z(v.z), w(v.w) {}

  Vec4<t>& operator=(const Vec4<t>& v) {
    if (this != &v) {
      x = v.x;
      y = v.y;
      z = v.z;
      w = v.w;
    }
    return *this;
  }

  Vec4<t> operator+(const Vec4<t>& v) const {
    return Vec4<t>(x + v.x, y + v.y, z + v.z, w + v.w);
  }
  Vec4<t> operator-(const Vec4<t>& v) const {
    return Vec4<t>(x - v.x, y - v.y, z - v.z, w - v.w);
  }
  Vec4<t> operator*(float f) const {
    return Vec4<t>(x * f, y * f, z * f, w * f);
  }
  t operator*(const Vec4<t>& v) const {
    return x * v.x + y * v.y + z * v.z + w * v.w;
  }

  float norm() const { return std::sqrt(x * x + y * y + z * z + w * w); }
  Vec4<t>& normalize(t l = 1) {
    *this = (*this) * (l / norm());
    return *this;
  }

  t& operator[](const int i) {
    if (i == 0) return x;
    if (i == 1) return y;
    if (i == 2) return z;
    if (i == 3) return w;
    throw std::out_of_range("Vec4 index out of range");
  }

  const t& operator[](const int i) const {
    if (i == 0) return x;
    if (i == 1) return y;
    if (i == 2) return z;
    if (i == 3) return w;
    throw std::out_of_range("Vec4 index out of range");
  }

  template <class>
  friend std::ostream& operator<<(std::ostream& s, Vec4<t>& v);
};

typedef Vec2<float> Vec2f;
typedef Vec2<int>   Vec2i;
typedef Vec3<float> Vec3f;
typedef Vec3<int>   Vec3i;
typedef Vec4<float> Vec4f;

template <> template <> Vec3<int>::Vec3(const Vec3<float>& v);
template <> template <> Vec3<float>::Vec3(const Vec3<int>& v);


template <class t> std::ostream& operator<<(std::ostream& s, Vec2<t>& v) {
    s << "(" << v.x << ", " << v.y << ")\n";
    return s;
}

template <class t> std::ostream& operator<<(std::ostream& s, Vec3<t>& v) {
    s << "(" << v.x << ", " << v.y << ", " << v.z << ")\n";
    return s;
}

//////////////////////////////////////////////////////////////////////////////////////////////

const int DEFAULT_ALLOC=4;

class Matrix {
    std::vector<std::vector<float> > m;
    int rows, cols;
public:
    Matrix(int r=DEFAULT_ALLOC, int c=DEFAULT_ALLOC);
    inline int nrows();
    inline int ncols();

    static Matrix identity(int dimensions);
    std::vector<float>& operator[](const int i);
    Matrix operator*(const Matrix& a);
    Matrix transpose();
    Matrix inverse();

    friend std::ostream& operator<<(std::ostream& s, Matrix& m);

    Vec4f operator*(const Vec4f& v) const {
      assert(rows == 4 && cols == 4);  // È·±£¾ØÕóÊÇ 4x4

      Vec4f result;
      result.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3] * v.w;
      result.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3] * v.w;
      result.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3] * v.w;
      result.w = m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3] * v.w;

      return result;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////


#endif __GEOMETRY_H__
