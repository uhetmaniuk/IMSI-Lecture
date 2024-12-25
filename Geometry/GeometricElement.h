#pragma once

#include <array>

#include "Utilities/Span.h"

namespace msfem {

enum GeometryType : char {Bar2 = 0, Tri3, Quad4, Tet4, Hex8};

/// \brief An abstract class to define a geometric element in a mesh.
struct GeometricElement {
    GeometricElement() = delete;
    
    virtual ~GeometricElement() = default;
    
    virtual GeometryType GetType() const = 0;
    
    virtual int Nodes() const = 0;
    
    virtual Span<int> NodesList() const = 0;
    
protected:

    GeometricElement(GeometryType type) : d_type(type) {}
    
    GeometryType d_type;

};

struct Bar: public GeometricElement {
    std::array<int, 2> d_nodes;

    Bar(int a, int b) : GeometricElement(GeometryType::Bar2),
    d_nodes({a, b}) {}
    
    ~Bar() = default;
    GeometryType GetType() const override { return GeometryType::Bar2; }
    int Nodes() const override { return 2; }
    Span<int> NodesList() const override { return {&array[0], 2}; }
};

struct Triangle : public GeometricElement {
    std::array<int, 3> d_nodes;

    Triangle(int a, int b, int c) : GeometricElement(GeometryType::Tri3),
    d_nodes({a, b, c}) {}
    
    ~Triangle() = default;
    GeometryType GetType() const override { return GeometryType::Tri3; }
    int Nodes() const override { return 3; }
    Span<int> NodesList() const override { return {&array[0], 3}; }
};

struct Quadrilateral : public GeometricElement {
    std::array<int, 4> d_nodes;
    
    Quadrilateral(int a, int b, int c, int d) : GeometricElement(GeometryType::Quad4),
    d_nodes({a, b, c, d}) {}
    
    ~Quadrilateral() = default;
    GeometryType GetType() const override { return GeometryType::Quad4; }
    int Nodes() const override { return 4; }
    Span<int> NodesList() const override { return {&array[0], 4}; }
};

struct Tetrahedron : public GeometricElement {
    std::array<int, 4> d_nodes;
    
    Tetrahedron(int a, int b, int c, int d) : GeometricElement(GeometryType::Tet4),
    d_nodes({a, b, c, d}) {}
    
    ~Tetrahedron() = default;
    GeometryType GetType() const override { return GeometryType::Tet4; }
    int Nodes() const override { return 4; }
    Span<int> NodesList() const override { return {&array[0], 4}; }
};

struct Hexahedron : public GeometricElement {
    std::array<int, 8> d_nodes;
    
    Hexahedron(int a, int b, int c, int d, int e, int f, int g, int h) : GeometricElement(GeometryType::Hex8),
    d_nodes({a, b, c, d, e, f, g, h}) {}
    
    ~Hexahedron() = default;
    GeometryType GetType() const override { return GeometryType::Hex8; }
    int Nodes() const override { return 8; }
    Span<int> NodesList() const override { return {&array[0], 8}; }
};

} // namespace msfem
