#pragma once
#pragma comment(lib, "libscip.dll")
#pragma comment(lib, "C:\\Program Files\\SCIPOptSuite 7.0.3\\lib\\libscip.lib")

#ifdef MATHLIBRARY_EXPORTS
#define MATHLIBRARY_API __declspec(dllexport)
#else
#define MATHLIBRARY_API __declspec(dllimport)
#endif

extern "C" MATHLIBRARY_API int write_features(const char* problem_filename, 
											  const char* csv_filename,
											  const char* instance_name);


