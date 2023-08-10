#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Search_traits_3.h>
#include <list>
#include <vector>
#include <iostream>
#include <sstream>
#include <ctime>
#include <stdio.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Search_traits_3<Kernel> TreeTraits;
typedef CGAL::Fair<TreeTraits> Fair;
typedef CGAL::Euclidean_distance<TreeTraits> Distance;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits, Distance, Fair> Neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

constexpr unsigned int JET_NEIGHBOR = 24;
constexpr unsigned int JET_NUM_FRAME = 40000;

inline void read_xyz_file(const std::string input_filename, std::list<Point>* points) {
	std::fstream in_file;
	in_file.open(input_filename, std::ios::in);
	std::string line, coord;

	while (getline(in_file, line))
	{
		std::stringstream s_line(line);
		float value[3]{};
		for (int i = 0; i != 3; ++i)
		{
			getline(s_line, coord, ' ');
			value[i] = std::stof(coord);
		}
		points->push_back(Point(value[0], value[1], value[2]));
	}
}

inline void write_xyz_file(const std::string output_filename, std::list<Point>* points) {
	std::fstream out_file;
	out_file.open(output_filename, std::ios::app);
	for (auto it = points->begin(); it != points->end(); ++it)
	{
		out_file << std::to_string(it->x()) << " "
			<< std::to_string(it->y()) << " "
			<< std::to_string(it->z()) << std::endl;
	}
	out_file.close();
}

int main() {
	/*const std::string input_filename = CGAL::data_file_path("D:\\Games\\data\\jet_test\\modified\\noisy\\0001.xyz");
	const std::string output_filename = "D:\\Games\\data\\jet_test\\modified\\jet\\0001.xyz";*/

	/*const std::string input_filename = CGAL::data_file_path("D:\\Games\\data\\bunny\\bun_add_noise.ply");
	const std::string output_filename = "D:\\Games\\data\\bunny\\bun_denoise.xyz";*/

	const std::string input_filename = CGAL::data_file_path("D:\\Games\\data\\TEST\\maps\\noise_0.05\\map_01.xyz");
	const std::string output_filename = "D:\\Games\\data\\TEST\\maps\\noise_0.05\\denoise_map_01.xyz";

	std::fstream f(output_filename.c_str());
	if (f.good())
	{
		std::cerr << "Error: File already exists!\n";
		return -1;
	}

	std::list<Point> points;

	if (!CGAL::IO::read_points(input_filename, std::back_inserter(points)))
	{
		std::cerr << "cant read\n";
		return -1;
	}

	const unsigned int points_size = (unsigned int)points.size();
	const unsigned int bucket_size = 2048;
	Fair fair(bucket_size);

	Tree tree(points.begin(), points.end(), fair);
	tree.statistics(std::cout);

	Neighbor_search search(tree, Point(0,0,0), points_size);
	std::cout << "sorted\n";

	unsigned int count(0);
	Neighbor_search::iterator start = search.begin();
	while (count < points_size - JET_NUM_FRAME) 
	{
		points.clear();
		for (Neighbor_search::iterator it = start; it != start + JET_NUM_FRAME; ++it)
		{
			points.push_back(it->first);
		}
		CGAL::jet_smooth_point_set<Concurrency_tag>(points, JET_NEIGHBOR);
		std::cout << "finish jet fitting with count " << count << std::endl;
		write_xyz_file(output_filename, &points);
		count += JET_NUM_FRAME;
		start += JET_NUM_FRAME;
	}

	points.clear();
	for (Neighbor_search::iterator it = start; it != search.end(); ++it)
	{
		points.push_back(it->first);
	}
	CGAL::jet_smooth_point_set<Concurrency_tag>(points, JET_NEIGHBOR);
	write_xyz_file(output_filename, &points);

	std::cout << "success\n";

	return 0;
}