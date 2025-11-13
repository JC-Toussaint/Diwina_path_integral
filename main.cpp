#include <iostream>
#include <signal.h>
#include <stdio.h>      // for perror()
#include <sys/stat.h>   // for mkdir(), stat()
#include <sys/types.h>  // for mkdir(), stat()
#include <unistd.h>     // for getpid(), stat()
#include <random>

#include "chronometer.h"
#include "fem.h"
#include "fem2d.h"
#include "mesh.h"

#include "pathIntegral.h"

#include <gmsh.h>
#include <execution>

inline double DegreesToRadians(double degrees) { return degrees * (M_PI / 180.0); }

Eigen::Vector3d rotation(Eigen::Vector3d M, Eigen::Vector3d axe, double theta_rad) {
	return Eigen::AngleAxis(theta_rad,axe) * M;
}

// Catch some deadly signals in order to save the state before quitting.
volatile sig_atomic_t received_signal = 0;

static void signal_handler(int signal_number) { received_signal = signal_number; }

// Create the output directory if it does not exist yet.
static void create_dir_if_needed(std::string dirname)
{
	const char *name = dirname.c_str();
	struct stat statbuf;
	int res = stat(name, &statbuf);
	if (res != 0 && errno != ENOENT)
	{
		std::cout << "could not be searched.\n";
		perror(name);
		exit(1);
	}
	if (res == 0)
	{  // path exists
		if (S_ISDIR(statbuf.st_mode))
		{
			std::cout << "(already exists)\n";
			return;
		}
		else
		{
			std::cout << "exists and is not a directory.\n";
			exit(1);
		}
	}

	// The directory does not exist (stat() reported ENOENT), create it.
	res = mkdir(name, 0777);
	if (res != 0)
	{
		std::cout << "could not be created.\n";
		perror(name);
		exit(1);
	}
	std::cout << "(created)\n";
}

// Return the number of characters in an UTF-8-encoded string.
static int char_length(const std::string &s)
{
	int len = 0;
	for (const char *p = s.c_str(); *p; ++p)
	{
		if ((*p & 0xc0) != 0x80)  // count bytes not matching 0x10xxxxxx
			len++;
	}
	return len;
}

// Return a string padded with spaces to a given length.
static std::string pad(const std::string &s, int length)
{
	length -= char_length(s);
	if (length < 0) length = 0;
	return s + std::string(length, ' ');
}

void prompt(void)
{
	std::cout << "\t┌────────────────────────────────┐\n";
	std::cout << "\t│          pathIntegral          │\n";
	std::cout << "\t│      CNRS Grenoble – INP       │\n";
	std::cout << "\t└────────────────────────────────┘\n";
}

std::string parseOptions(Settings &settings, int argc, char *argv[], unsigned int &random_seed)
{
	int print_help = false;
	int print_version = false;
	int print_defaults = false;
	int verify = false;
	int use_fixed_seed = false;

	struct Option
	{
		std::string short_opt, long_opt;
		const char *help;
		int *setting;
	};
	struct Option options[] = {
			{"-h", "--help", "display short help and exit", &print_help},
			{"-V", "--version", "display version information and exit", &print_version},
			{"", "--print-defaults", "print default settings and exit", &print_defaults},
			{"", "--verify", "verify a settings file and exit", &verify},
			{"-v", "--verbose", "enable verbose mode", &settings.verbose},
			{"", "--seed", "set random seed", &use_fixed_seed},
			{"", "", nullptr, nullptr}  // sentinel
	};

	int optind;
	for (optind = 1; optind < argc; optind++)
	{
		char *opt = argv[optind];
		Option *o;
		for (o = options; o->setting; o++)
		{
			if (opt == o->short_opt || opt == o->long_opt)
			{
				(*o->setting)++;
				if (o->long_opt == "--seed" && optind < argc - 1)
				{
					random_seed = atol(argv[++optind]);
				}
				break;
			}
		}
		if (!o->setting)  // option not found
			break;
	}
	if (print_help)
	{
		std::cout << "Usage: " << argv[0] << " [options] settings_file\n";
		std::cout << "Options:\n";
		for (Option *o = options; o->setting; o++)
		{
			std::cout << "  ";
			if (!o->short_opt.empty())
				std::cout << o->short_opt << " ";
			else
				std::cout << "   ";
			std::cout << pad(o->long_opt, 17) << " " << o->help << "\n";
		}
		exit(0);
	}
	if (print_version)
	{
		std::cout << "pathIntegral " << pathIntegral_version << "\n";
		exit(0);
	}
	if (print_defaults)
	{
		Settings::dumpDefaults();
		exit(0);
	}
	if (optind != argc - 1)
	{
		std::cerr << "Usage: pathIntegral [options] config_file.yml\n";
		exit(1);
	}
	std::string filename = argv[optind];
	if (verify)
	{
		settings.read(filename);
		settings.infos();
		exit(0);
	}
	if (!use_fixed_seed)
	{
		// Use a truly random seed.
		std::random_device rd;
		random_seed = rd();
	}

	return filename;
}

int main(int argc, char *argv[])
{
	Settings mySettings;
	chronometer counter;

	unsigned int random_seed;
	std::string filename = parseOptions(mySettings, argc, argv, random_seed);
	srand(random_seed);
	prompt();
	
	if (mySettings.verbose) std::cout << "verbose mode:      on\n";
	std::cout << "pathIntegral version: " << pathIntegral_version << '\n';
	if (ENABLE_SEQ) std::cout << "sequential mode:   on\n";
	else std::cout << "parallel mode:     on\n";
	std::cout << "process ID:        " << std::to_string(getpid()) << '\n';
	std::cout << "random seed:       " << random_seed << '\n';
	mySettings.setFileDisplayName(filename == "-" ? "standard input" : filename);
	std::cout << "settings file:     " << mySettings.getFileDisplayName() << '\n';
	if (!mySettings.read(filename))
	{
		std::cerr << "Error: no settings found.\n";
		return 1;
	}
	std::cout << "mesh file:         " << mySettings.getPbName() << '\n';
	std::cout << "output directory:  " << mySettings.r_path_output_dir << " ";
	create_dir_if_needed(mySettings.r_path_output_dir);
	Fem fem = Fem(mySettings);

	if (mySettings.verbose)
	{
		std::cout << "-- settings: -----------------------------------\n";
		mySettings.infos();
		std::cout << "-- end of settings -----------------------------\n";
		fem.msh.infos();
	}

	counter.reset();
	std::cout << "starting on:       " << date() << std::endl;

	chronometer fmm_counter(2);

	// Catch SIGINT and SIGTERM.
	struct sigaction action;
	action.sa_handler = signal_handler;
	action.sa_flags = 0;
	sigemptyset(&action.sa_mask);
	if (sigaction(SIGINT, &action, NULL) == -1)
	{
		perror("SIGINT");
		return EXIT_FAILURE;
	}
	if (sigaction(SIGTERM, &action, NULL) == -1)
	{
		perror("SIGTERM");
		return EXIT_FAILURE;
	}

    mySettings.p_rot.infos();

	Eigen::Vector3d axe1{mySettings.p_rot.axe1.data()};
	double angle1 =  DegreesToRadians(mySettings.p_rot.angle1);

	Eigen::Vector3d axe2{mySettings.p_rot.axe2.data()};
	double angle2 =  DegreesToRadians(mySettings.p_rot.angle2);

	for (int i=0; i<fem.msh.getNbNodes(); i++)
	    {
		Eigen::Vector3d xyz= fem.msh.getNode_p(i);
		xyz = rotation(xyz, axe1, angle1);
		xyz = rotation(xyz, axe2, angle2);
		fem.msh.setNode_p(i, xyz);

		Eigen::Vector3d u= fem.msh.getNode_u(i);
		u = rotation(u, axe1, angle1);
		u = rotation(u, axe2, angle2);
		fem.msh.setNode_u(i, u);
	    }

	fem.msh.update_lengths();
	
	std::cout << "sizes  : " << fem.msh.l << std::endl;
	std::cout << "filled : " << (mySettings.filled ? "true" : "false") << std::endl;

    PathInt::pathIntegral myStuff(fem.msh,mySettings.filled);

	Fem2d fem2d(fem.msh.c, fem.msh.l, mySettings.p_electrostatics.CE,   mySettings.p_electrostatics.V, 
					  mySettings.p_detector.zoomFactor, mySettings.p_detector.meshSize);
	fem2d.infos();

	std::vector<int> nodeIndices(fem2d.getNbNodes());
	std::iota(nodeIndices.begin(), nodeIndices.end(), 0);

	std::for_each(std::execution::par, nodeIndices.begin(), nodeIndices.end(),
	              [&](int nod) { myStuff.processNode(nod, fem2d, mySettings); });

	// exporting Magnetization integrals along the beam
 	fem2d.exportMagIntegrals(mySettings.getSimName());

    fem2d.exportRatioGrayScaleImage(mySettings, ExportType::CONTRAST);
	fem2d.exportRatioGrayScaleImage(mySettings, ExportType::MZ_INTEGRAL);
	fem2d.exportRatioGrayScaleImage(mySettings, ExportType::PATH_LENGTH);

	double total_time = counter.fp_elapsed();
	std::cout << "\nComputing time: " << counter.convertSeconds(total_time);

        for (auto& node : fem2d.getNodes()){ //recuperation des numeros et aimantations integrees pour chaque noeud du maillage 2d
                /* M' = ez x M  => Mx'=-My  My'=Mx */
                node.Mx = -node.My_integral;
                node.My = +node.Mx_integral;
        }
	fem2d.util();

	std::cout << "Fast Multipole Calculation\n";
	// int ierr = pot2D::fmm2d_sum(fem2d);
    //int ierr = pot2D::scalfmm2d_sum(fem2d);
    int ierr = pot2D::direct2d_sum(fem2d);
    std::cout << "fmm2D returned " << ierr << std::endl;

 	fem2d.exportHoloPhase(mySettings.getSimName());
 	fem2d.exportRatioRGBscaleImage(mySettings, ExportType::HOLO_PHASE);
	return 0;
}
