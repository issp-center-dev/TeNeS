if [ -z ${BOOST_ROOT} ]; then
  BOOST_ROOT=/opt/homebrew/include
fi

rm -rf include
mkdir include

bcp --boost=${BOOST_ROOT} \
  core/demangle.hpp \
  optional \
  multiprecision/cpp_int.hpp \
  \
  include
