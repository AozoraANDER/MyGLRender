#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include "geometry.h"
#include "model.h"
#include "tgaimage.h"

#define M_PI 3.14159265358979323846

// 计算叉积函数
Vec3f cross(const Vec3f &v1, const Vec3f &v2) {
  return Vec3f(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
               v1.x * v2.y - v1.y * v2.x);
}

// 计算重心坐标
Vec3f barycentric(Vec2i A, Vec2i B, Vec2i C, Vec2i P) {
  Vec3f s[2];
  s[0] = Vec3f(C.x - A.x, B.x - A.x, A.x - P.x);
  s[1] = Vec3f(C.y - A.y, B.y - A.y, A.y - P.y);
  Vec3f u = cross(s[0], s[1]);  // 使用叉积函数计算

  if (std::abs(u.z) < 1) return Vec3f(-1, 1, 1);  // 三角形退化的情况

  return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
}

// 绘制三角形，使用深度缓冲
void triangle(Vec3f *world_coords, Vec2i *screen_coords, float *zbuffer,
              TGAImage &image, TGAColor color) {
  Vec2f bboxmin(std::numeric_limits<float>::max(),
                std::numeric_limits<float>::max());
  Vec2f bboxmax(-std::numeric_limits<float>::max(),
                -std::numeric_limits<float>::max());
  Vec2f clamp(image.get_width() - 1, image.get_height() - 1);

  // 计算三角形的包围盒
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      bboxmin[j] =
          std::max(0.f, std::min(bboxmin[j], (float)screen_coords[i][j]));
      bboxmax[j] =
          std::min(clamp[j], std::max(bboxmax[j], (float)screen_coords[i][j]));
    }
  }

  Vec3f P;
  for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
    for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
      Vec3f bc_screen = barycentric(screen_coords[0], screen_coords[1],
                                    screen_coords[2], Vec2i(P.x, P.y));
      if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;

      // 计算 P 点的深度值 (Z 值)
      P.z = 0;
      for (int i = 0; i < 3; i++) {
        P.z += world_coords[i].z * bc_screen[i];
      }

      int idx = int(P.x) + int(P.y) * image.get_width();
      if (zbuffer[idx] < P.z) {
        zbuffer[idx] = P.z;
        image.set(P.x, P.y, color);
      }
    }
  }
}

int main() {
  TGAImage image(800, 800, TGAImage::RGB);
  int width = 800;
  int height = 800;
  Model model("assets/african_head.obj");

  float *zbuffer = new float[width * height];
  for (int i = 0; i < width * height; ++i) {
    zbuffer[i] = std::numeric_limits<float>::lowest();  // 初始设为负无穷
  }

  Vec3f light_dir(0, 0, -1);  // 定义光照方向

  for (int i = 0; i < model.nfaces(); i++) {
    std::vector<int> face = model.face(i);
    Vec3f world_coords[3];
    Vec2i screen_coords[3];  // 使用 Vec2i，因为屏幕坐标是二维的

    for (int j = 0; j < 3; j++) {
      Vec3f v = model.vert(face[j]);
      screen_coords[j] =
          Vec2i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2.);
      world_coords[j] = v;  // 保存 3D 世界坐标
    }

    Vec3f n = cross(world_coords[2] - world_coords[0],
                    world_coords[1] - world_coords[0]);
    n.normalize();

    float intensity = n * light_dir;  // 计算光照强度

    if (intensity > 0) {
      triangle(
          world_coords, screen_coords, zbuffer, image,
          TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
    }
  }

  image.flip_vertically();  // 翻转图像，使原点在左下角
  image.write_tga_file("Alpha.tga");

  delete[] zbuffer;  // 释放 zbuffer
  return 0;
}
