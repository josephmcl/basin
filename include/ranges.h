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
            index(index), Linrange(Linrange) {}

        T &operator *() { return this->value(); }
        std::add_const_t<T &> operator * () const { this->value(); } 
        I &operator ++ () { /* weakly_incrementable */ 
            _inc_impl(); return *this; } 
        I operator ++(int) { /* incrementable */ 
            I temp = *this; ++*this; return temp; }     
        I &operator --() { /* weakly_incrementable */ 
            _dec_impl(); return *this; }             
        I operator --(int) { /* incrementable */ 
            I temp = *this; --*this; return temp; }  
        I &operator = (I const &lhs); /* assignable_from */
        bool operator == (const I &other) /* equality_comparable */
            const { return index == other.index; }
        bool operator != (const I &other) const { 
            return !(*this == other); } 

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
            std::size_t(_size)), this}; };

    
};  

    template <typename T, typename R>
    auto operator |(linrange<T> t, R r) {
        t.transforms.push_back(r);
        return t;
    }

}
