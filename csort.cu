#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/scan.h>

extern "C" {
  void sort_int_wrapper( int *data, int N)
  {
    thrust::device_ptr<int> dev_ptr(data);
    thrust::sort(dev_ptr, dev_ptr+N);
  }
  void sort_double_wrapper( double *data, int N)
  {
    thrust::device_ptr<double> dev_ptr(data);
    thrust::sort(dev_ptr, dev_ptr+N);
  }
  void sort_by_key_int_wrapper( int *key, int N, int *data)
  {
    thrust::device_ptr<int> dev_key(key);
    thrust::device_ptr<int> dev_ptr(data);
    thrust::stable_sort_by_key(dev_key, dev_key+N, dev_ptr);
  }
  void sort_by_key_double_wrapper( double *key, int N, int *data)
  {
    thrust::device_ptr<double> dev_key(key);
    thrust::device_ptr<int> dev_ptr(data);
    thrust::stable_sort_by_key(dev_key, dev_key+N, dev_ptr);
  }
  void prefixSum_wrapper( int *data, int *out ,int N)
  {
    thrust::device_ptr<int> dev_data(data);
    thrust::device_ptr<int> dev_out(out);
    thrust::exclusive_scan(dev_data, dev_data+N, dev_out, 1);
  }
}
