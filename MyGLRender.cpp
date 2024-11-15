#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include <vector> 
#include <algorithm>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <limits>

#define M_PI 3.14159265358979323846

std::vector<Vec2f> uv_coords;

Matrix perspective_projection(float fov, float aspect_ratio, float near,
                              float far) {
  float f = 1.0f / std::tan(fov / 2.0f);
  Matrix projection = Matrix::identity(4);

  projection[0][0] = f / aspect_ratio;
  projection[1][1] = f;
  projection[2][2] = (far + near) / (near - far);
  projection[2][3] = (2 * far * near) / (near - far);
  projection[3][2] = -1;
  projection[3][3] = 0;

  return projection;
}

Vec3f cross(const Vec3f &v1, const Vec3f &v2) {
  return Vec3f(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
               v1.x * v2.y - v1.y * v2.x);
}

// 计算叉积函数
void loadUVCoords(const std::string &filename, std::vector<Vec2f> &uv_coords) {
  std::ifstream in(filename, std::ifstream::in);
  if (in.fail()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream iss(line);
    char trash;
    if (!line.compare(0, 3, "vt ")) {
      iss >> trash >> trash;  // skip "vt"
      Vec2f uv;
      iss >> uv.x >> uv.y;
      uv_coords.push_back(uv);
      std::cout << "Loaded UV coordinate: (" << uv.x << ", " << uv.y << ")"
                << std::endl;
    }
  }
}

std::vector<std::vector<int>> face_uv_indices;

void loadFaceUVIndices(const std::string &filename,
                       std::vector<std::vector<int>> &face_uv_indices) {
  std::ifstream in(filename, std::ifstream::in);
  if (in.fail()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream iss(line);
    char trash;
    if (!line.compare(0, 2, "f ")) {
      iss >> trash;  // skip 'f'
      std::vector<int> face_uv;
      int v_idx, t_idx, n_idx;
      for (int i = 0; i < 3; i++) {
        char separator;
        iss >> v_idx >> separator >> t_idx >> separator >> n_idx;
        face_uv.push_back(t_idx - 1);  // 索引从1开始，改为从0开始
      }
      face_uv_indices.push_back(face_uv);
      std::cout << "Face UV indices: " << t_idx << std::endl;
    }
  }
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
void triangle(Vec3f *world_coords, Vec2i *screen_coords, Vec2f *uv_coords,
              float *zbuffer, TGAImage &image, TGAImage &diffuse_map,
              float intensity) {
  // 包围盒计算和深度缓冲代码保持不变

    Vec2f bboxmin(std::numeric_limits<float>::max(),
                std::numeric_limits<float>::max());
  Vec2f bboxmax(-std::numeric_limits<float>::max(),
                -std::numeric_limits<float>::max());
  Vec2f clamp(image.get_width() - 1, image.get_height() - 1);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      bboxmin[j] = std::max(
          0.f, std::min(bboxmin[j], static_cast<float>(screen_coords[i][j])));
      bboxmax[j] = std::min(
          clamp[j],
          std::max(bboxmax[j], static_cast<float>(screen_coords[i][j])));

    }
  }

  // 在三角形内采样纹理颜色
  

  // 计算 P 点的深度值 (Z 值)
  Vec3f P;
  for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
    for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
      // 使用barycentric重心坐标的部分
      Vec3f bc_screen = barycentric(screen_coords[0], screen_coords[1],
                                    screen_coords[2], Vec2i(P.x, P.y));
      if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;

      // 计算 P 点的深度值 (Z 值)
      P.z = 0;
      for (int i = 0; i < 3; i++) {
        P.z += world_coords[i].z * bc_screen[i];
      }

      // 更新zbuffer并设置像素颜色
      int idx = static_cast<int>(P.x) + static_cast<int>(P.y) * image.get_width();
      if (zbuffer[idx] < P.z) {
        zbuffer[idx] = P.z;

        // 计算纹理坐标
        Vec2f uv(0, 0);
        for (int i = 0; i < 3; i++) {
          uv.x += uv_coords[i].x * bc_screen[i];
          uv.y += uv_coords[i].y * bc_screen[i];
        }
        uv.x *= diffuse_map.get_width();
        uv.y *= diffuse_map.get_height();

        // 应用光照强度
        TGAColor color =
            diffuse_map.get(static_cast<int>(uv.x), static_cast<int>(uv.y));
        color.b = static_cast<unsigned char>(color.b * intensity);
        color.g = static_cast<unsigned char>(color.g * intensity);
        color.r = static_cast<unsigned char>(color.r * intensity);

       Vec3f bc_screen = barycentric(screen_coords[0], screen_coords[1],
                                      screen_coords[2], Vec2i(P.x, P.y));

        if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) {
         return;  // 退出当前函数
         }
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

  loadUVCoords("assets/african_head.obj", uv_coords);
  loadFaceUVIndices("assets/african_head.obj", face_uv_indices);

  TGAImage diffuse_map;
  if (!diffuse_map.read_tga_file("assets/african_head_diffuse.tga")) {
    std::cerr << "Error: Cannot load the texture file!" << std::endl;
    return 1;
  }
  diffuse_map.flip_vertically();

  float *zbuffer = new float[width * height];
  for (int i = 0; i < width * height; ++i) {
    zbuffer[i] = std::numeric_limits<float>::lowest();  // 初始设为负无穷
  }

  Vec3f light_dir(0, 0, -1);  // 定义光照方向

  // 定义透视投影矩阵
  float fov = 45* M_PI / 180.0f;
  float aspect_ratio = width / (float)height;
  float near = 0.1f;
  float far = 100.0f;
  Matrix projection = perspective_projection(fov, aspect_ratio, near, far);

  for (int i = 0; i < model.nfaces(); i++) {
    std::vector<int> face = model.face(i);
    std::vector<int> face_uv = face_uv_indices[i];  // 获取每个面的UV索引

    Vec3f world_coords[3];
    Vec2i screen_coords[3];
    Vec2f uv_coords_triangle[3];

    for (int j = 0; j < 3; j++) {
      Vec3f v = model.vert(face[j]);

      // 将 Vec3f 扩展为 Vec4f，用于与投影矩阵相乘

      Vec4f homogenous_vertex(v.x, v.y, v.z + 3.0f, 1.0f);  // 将模型沿 Z 轴后移

      Vec4f transformed_vertex = projection * homogenous_vertex;



      // 进行透视除法
      if (transformed_vertex.w != 0) {
        transformed_vertex.x /= transformed_vertex.w;
        transformed_vertex.y /= transformed_vertex.w;
        transformed_vertex.z /= transformed_vertex.w;
      }

      // 转换为屏幕坐标
      screen_coords[j] = Vec2i((transformed_vertex.x + 1.0f) * width / 2.0f,
                               (transformed_vertex.y + 1.0f) * height / 2.0f);

      world_coords[j] = v;
      uv_coords_triangle[j] = uv_coords[face_uv[j]];  // 获取每个顶点的UV坐标
    }


    Vec3f n = cross(world_coords[2] - world_coords[0],
                    world_coords[1] - world_coords[0]);

    n.normalize();
    float intensity = n * light_dir;  // 计算光照强度

    if (intensity > 0) {
      triangle(world_coords, screen_coords, uv_coords_triangle, zbuffer, image,
               diffuse_map, intensity);

    }
  }

  //image.flip_vertically();  // 翻转图像，使原点在左下角
  image.write_tga_file("Alpha.tga");

  delete[] zbuffer;  // 释放 zbuffer
  return 0;
}
