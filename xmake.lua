set_xmakever("3.0.0")

--add_requires("openmp", { optional = true })
add_rules("mode.debug", "mode.release")

local src_files = {
	"src/main.cxx",
	"src/io/*.cxx",
	"src/func/*.cxx",
	"src/basic/*.cxx",
}

target("pkcpp")
set_kind("binary")
add_files(src_files)

add_includedirs("ext/tomlpp/include", { public = false })

set_languages("c++23")

set_targetdir("build")

set_policy("build.c++.modules", true)

add_packages("openmp")

--on_load(function(target)
--	if target:pkg("openmp") then
--		cprint("${green}OpenMP found: build in parallel mode.")
--	else
--		cprint("${yellow}OpenMP not found: Build in serial mode.")
--	end
--end)
