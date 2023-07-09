# https://gitee.com/haha-web/eigen.git

# 引入doctest
include(FetchContent)
FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitee.com/haha-web/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(Eigen)