Developer guide
===============

Here we have a few conventions for developing in `libcanon`.

1. All code are written in *modern* C++. [C++ core guidelines][Core] are
   generally followed.

2. C++14 compatibility is generally emphasized, although the concept-lite TS is
   used in the code in an optional way.

3. The coding style basically follows the [WebKit code style
   guidelines][WebKit], with a few exceptions,

   1. Naming basically follows the convention in Stroutrup [PPP style
      guide][PPP].  

   2. Member names and accessors are named according to the [Google C++ style
      guide][Google].

4. For functions inputs, generally references are used when the function does
   not do anything with the deallocation of the object, pointers can be used
   when no value, or `NULL` is a possibility, and unique pointers are used when
   the function takes the ownership of the pointer.  Similar principles applies
   for function return values.

5. Parenthesis initialization is preferred over brace initialization, which is
   only used when it is necessary.

6. The code base should be compatible with both G++ later than 6.3.0 and
   clang++ later than 3.8.

[Core]: https://github.com/isocpp/CppCoreGuidelines
[WebKit]: https://webkit.org/code-style-guidelines/
[PPP]: http://www.stroustrup.com/Programming/PPP-style-rev3.pdf
[Google]: https://google.github.io/styleguide/cppguide.html

