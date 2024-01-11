#pragma once 
/*  range.h  
    
    Implements a range object similar to the range functions found in 
    many languages and libraries. 

    Initial design take from stackoverflow.com/a/56217827 */

#include <ranges>
#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>

namespace numerics {

/*  template <typename T> struct numerics::range 

    For a given range {from, to} and a regular spacing, generate an
    iterator that iterates from "from" to "to" on a step length of 
    "space." */ 
template <typename T> struct range {
    T from, space, to;
    range(T from, T space, T to): from(from), space(space), to(to) {}
    struct iterator {
        const T to, space;
        T current;
        T operator *() { return current; } 
        iterator &operator ++() { 
            /* Someone will probably say this is a no-no. Oh well. */
            current += space; 
            if(current >= to)
                current = to;
            return *this; }
        bool operator ==(const iterator &other) const { 
            return current == other.current; }
        bool operator !=(const iterator &other) const { 
            return !(*this == other); } };
    iterator begin() const { return iterator{to, space, from}; }
    iterator end()   const { return iterator{to, space,   to}; } };

/*  template <typename T> struct numerics::linrange 

    For a given range {from, to} with n nodes, generate an iterator 
    that yields an evenly spaced range (h = to - from / n - 1). */

template <typename T> struct linrange {
    std::size_t const _size;
    T const from, to, space;
    std::vector<T(*)(T)> transforms = {};
    linrange() = default;
    linrange(T from, T to, std::size_t size): _size(size), from(from), 
        to(to), space((to - from) / static_cast<T>(size - 1)) {}
    linrange(linrange<T> const &lr, bool transforms=true): 
        _size(lr._size), from(lr.from), to(lr.to), space(lr.space), 
        transforms(transforms? lr.transforms: std::vector<T(*)(T)>()){}
    struct iterator {
        
        std::size_t index;
        linrange const *Linrange;
        T _value;

        /* Several of the below definitions fulfill constraints for:
            
            - std::input_iterator 
            - std::forward_iterator */

        using I = iterator; 
        using value_type = T; /* indirectly_readable */                
        using difference_type = std::ptrdiff_t;

        iterator() {}; /* constructible_from */
        iterator(std::size_t index, linrange const *Linrange): 
            index(index), Linrange(Linrange) { value(); }

        T &operator *() { return this->value(); }
        std::add_const_t<T &> operator *() const { this->value(); } 
        I &operator ++() { /* weakly_incrementable */ 
            _inc_impl(); return *this; } 
        I operator ++(int) { /* incrementable */ 
            I temp = *this; ++*this; return temp; }     
        I &operator --() { /* weakly_incrementable */ 
            _dec_impl(); return *this; }             
        I operator --(int) { /* incrementable */ 
            I temp = *this; --*this; return temp; }  
        //inline I &operator = (Type* rhs) {_ptr = rhs; return *this;}
        inline I &operator = (I const &that) {index = that.index; return *this;}
        /* assignable_from */
        bool operator == (const I &that) /* equality_comparable */
            const { return index == that.index; }
        bool operator != (const I &that) const { 
            return !(index == that.index); } 

        bool operator < (I const &that) const { 
            return _value < that._value; }

        bool operator > (I const &that) const { 
            return _value > that._value; } 

        bool operator <= (I const &that) const { 
            return _value <= that._value; }

        bool operator >= (I const &that) const { 
            return _value >= that._value; } 

        inline std::iter_difference_t<I> operator - (I const &that) const { 
            return iterator{index - that.index, Linrange}; } 

        inline I operator - (std::size_t const &n) const { 
            return iterator{index - n, Linrange}; } 

        inline I &operator -= (std::size_t const &that) { 
            index -= that; return *this; } 

        friend inline I 
        operator - (const std::size_t &lhs, const I &rhs) {
            return iterator{lhs + rhs.index, rhs.Linrange}; } 

        inline std::iter_difference_t<I> 
        operator + (I const &that) const { 
            return iterator{index + that.index, Linrange}; } 

        inline I 
        operator + (std::size_t const &n) const { 
            return iterator{index + n, Linrange}; } 

        inline I &operator += (std::size_t const &that) { 
            index += that; return *this; } 

        friend inline I 
        operator + (const std::size_t &lhs, const I &rhs) {
            return iterator{lhs + rhs.index, rhs.Linrange}; } 

        inline std::iter_reference_t<I> 
        operator[](const std::size_t &n) const { 
            index = n; return *this;}

        void _inc_impl() { 
            index += 1; 
            if(index >= Linrange->_size)
                index = Linrange->_size; } 

        void _dec_impl() { 
            index -= 1; 
            if(index <= 0)
                index = 0; } 

        T &value() {
            _value = Linrange->from + (index * Linrange->space); 
            for (auto &t : Linrange->transforms)
                _value = t(_value);  
            return _value; }

    };
    
    /* std::forward_range constraints */
    auto begin() const noexcept { return iterator{0, this}; }
    auto end()   const noexcept { return iterator{_size, this}; }
    /* std::sized_range constraint */                             
    auto size()  const noexcept { return _size; }; 
    auto operator [](std::size_t const index) const noexcept {
        return iterator{std::clamp(index, std::size_t(0), 
            std::size_t(_size)), this}; 
    };

    template <typename R>
    linrange<T> &operator |= (R r) {
        transforms.push_back(r);
        return *this;
    }

    linrange<T> clip(std::size_t start, std::size_t stop) const {
        T to   = this->from + ((stop + 1)  * space); 
        T from = this->from + (start * space); 
        std::size_t size = _size - start - (_size - (stop + 1));
        auto that = linrange(from, to, size);
        that.transforms = transforms;
        return that;
    }

    iterator min_element() {
        auto first = begin();
        auto last  = end();
        if (first == last) return last;
        iterator smallest = first;
        ++first;
        for (; first != last; ++first) {
            if (*first < *smallest)
                    smallest = first;
        }
        return smallest;
    }
    
};  

    template <typename T, typename R>
    auto operator |(linrange<T> t, R r) {
        t.transforms.push_back(r);
        return t;
    }
}
