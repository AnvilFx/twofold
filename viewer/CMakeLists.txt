include(FetchContent)

set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_GAMES OFF CACHE BOOL "" FORCE)
set(CUSTOMIZE_BUILD OFF CACHE BOOL "" FORCE)
FetchContent_Declare(raylib
    URL https://github.com/raysan5/raylib/archive/refs/tags/5.5.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_Declare(imgui
    URL https://github.com/ocornut/imgui/archive/refs/tags/v1.92.7.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_Declare(rlimgui
    URL https://github.com/raylib-extras/rlImGui/archive/118221c8a532501a3c2c236b62950877cbefc2b4.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(raylib imgui rlimgui)

add_library(viewer_imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${rlimgui_SOURCE_DIR}/rlImGui.cpp
)
target_include_directories(viewer_imgui PUBLIC ${imgui_SOURCE_DIR} ${rlimgui_SOURCE_DIR})
target_link_libraries(viewer_imgui PUBLIC raylib)

add_executable(twofold_viewer draw.cpp main.cpp scenes.cpp session.cpp)
target_link_libraries(twofold_viewer PRIVATE twofold viewer_imgui)
if(WIN32)
    set_target_properties(twofold_viewer PROPERTIES WIN32_EXECUTABLE TRUE)
    if(MSVC)
        target_link_options(twofold_viewer PRIVATE /ENTRY:mainCRTStartup)
    endif()
endif()
