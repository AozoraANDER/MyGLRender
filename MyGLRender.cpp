#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include <vector> 
#include <iostream> 

#define M_PI 3.14159265358979323846

Vec3f barycentric(Vec2i A, Vec2i B, Vec2i C, Vec2i P) {
  Vec3f s[2];
  s[0] = Vec3f(C.x - A.x, B.x - A.x, A.x - P.x);
  s[1] = Vec3f(C.y - A.y, B.y - A.y, A.y - P.y);
  Vec3f u = s[0] ^ s[1];  // 计算叉积

 
  if (std::abs(u.z) < 1) return Vec3f(-1, 1, 1);

  return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
}

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
  // 找到边界框的边界
  int minX = std::min(t0.x, std::min(t1.x, t2.x));
  int minY = std::min(t0.y, std::min(t1.y, t2.y));
  int maxX = std::max(t0.x, std::max(t1.x, t2.x));
  int maxY = std::max(t0.y, std::max(t1.y, t2.y));

  // 遍历边界框中的每一个像素
  for (int y = minY; y <= maxY; y++) {
    for (int x = minX; x <= maxX; x++) {
      Vec2i P(x, y);
      Vec3f bary = barycentric(t0, t1, t2, P);

      // 如果重心坐标的所有分量都是非负，说明点在三角形内
      if (bary.x >= 0 && bary.y >= 0 && bary.z >= 0) {
        image.set(x, y, color);  // 填充颜色
      }
    }
  }
}

int main() {
  TGAImage image(800, 800, TGAImage::RGB);
  int width = 800;
  int height = 800;
  Model model("assets/african_head.obj");

  Vec3f light_dir(0, 0, -1);  // define light_dir

  for (int i = 0; i < model.nfaces(); i++) {
    std::vector<int> face = model.face(i);
    Vec2i screen_coords[3];
    Vec3f world_coords[3];
    for (int j = 0; j < 3; j++) {
      Vec3f v = model.vert(face[j]);
      screen_coords[j] =
          Vec2i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2.);
      world_coords[j] = v;
    }
    Vec3f n = (world_coords[2] - world_coords[0]) ^
              (world_coords[1] - world_coords[0]);
    n.normalize();
    float intensity = n * light_dir;
    if (intensity > 0) {
      triangle(
          screen_coords[0], screen_coords[1], screen_coords[2], image,
          TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
    }
  }
  image.flip_vertically();
  image.write_tga_file("Alpha.tga");
  return 0;
}