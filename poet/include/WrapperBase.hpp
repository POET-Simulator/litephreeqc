#pragma once

#include <phrqtype.h>
#include <span>

class WrapperBase {
public:
  virtual ~WrapperBase() = default;

  std::size_t size() const { return this->num_elements; };

  virtual void get(std::span<LDBLE> &data) const = 0;

  virtual void set(const std::span<LDBLE> &data) = 0;

protected:
  std::size_t num_elements = 0;
};