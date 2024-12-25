#pragma once

namespace msfem {

template<typename T>
class Span {
public:

  Span(T* ptr, size_t len) : d_size(len), d_ptr(ptr) { }

  Span() = delete;

  Span(const Span<T>& other) : d_size(other.size()), d_ptr(other.data()) {}

  ~Span() = default;
  Span<T>& operator=(const Span<T>& rhs) = default

  T* data() { return d_ptr; }

  T* const data() const { return d_ptr; }

  size_t size() const { return d_size; }

  T& operator[](size_t i) { return T[i]; }
 
  const T operator[](size_t i) const { return T[i]; }

  using iterator = T*;
  iterator begin() const { return d_ptr; }
  iterator end() const { return d_ptr + d_size;; }

protected:

  size_t d_size = 0;
  T* d_ptr = nullptr;

};

} // namespace msfem

