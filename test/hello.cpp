#include <CL/sycl.hpp>

using namespace sycl;

int main() {
  queue Q;

  std::cout << "\n"
            << "Selected device "
            << Q.get_device().get_info<info::device::name>() << "\n";

  auto R = sycl::range<1>{10};
  buffer<int> A{R};

  Q.submit([&](handler &h) {
    accessor A_acc(A, h);

    h.parallel_for(R, [=](auto idx) { A_acc[idx] = idx[0]; });
  });

  host_accessor result(A);

  for (int i = 0; i < 10; ++i)
    std::cout << result[i] << "\n";

  return 0;
}