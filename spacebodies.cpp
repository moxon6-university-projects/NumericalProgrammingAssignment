#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>

using namespace std;

struct Body {

    // Position of Body
    double x;
    double y;
    double z;

    //Velocity of Body
    double vx;
    double vy;
    double vz;

    //Acceleration of Body
    double ax;
    double ay;
    double az;

    //Mass of Body
    double mass;

    //Does body still 'exist' ( active=False => Took part in merging process earlier and is now ignored)
    bool active;
};

//Body Setup Functions
void setup_along_line(Body *bodies, int number_bodies) {
    for (int i = 0; i < number_bodies; i++) {
        struct Body new_body = {1000 * i, 0.0, 0.0,  //Position
                                0.0, 0.0, 0.0,  //Velocity
                                0.0, 0.0, 0.0,  //Acceleration
                                1.0,            //Mass
                                true};          //Initially Active
        bodies[i] = new_body;
    }
}

void setup_random_cube(Body *bodies, int number_bodies, double cube_width) {
    double x;
    double y;
    double z;

    for (int i = 0; i < number_bodies; i++) {

        x = cube_width / 2 + (cube_width / 2 * ((double) rand()) / (RAND_MAX));
        y = cube_width / 2 + (cube_width / 2 * ((double) rand()) / (RAND_MAX));
        z = cube_width / 2 + (cube_width / 2 * ((double) rand()) / (RAND_MAX));

        struct Body new_body = {x, y, z,  //Position
                                x / cube_width, y / cube_width, z / cube_width,  //Velocity
                                0.0, 0.0, 0.0,  //Acceleration
                                1.0,            //Mass
                                true};          //Initially Active
        bodies[i] = new_body;
    }
}

void setup_random_sphere(Body *bodies, int number_bodies, double sphere_radius) {
    struct Body new_body = {0, 0, 0,
                            0, 0, 0,
                            0, 0, 0,
                            1000000,
                            true};
    bodies[0] = new_body;
    double theta;
    double phi;
    for (int i = 1; i < number_bodies; i++) {

        //Random Angles
        theta = (double) 2 * M_PI * rand() / (RAND_MAX);
        phi = (double) 2 * M_PI * rand() / (RAND_MAX);

        //Polar -> Cartesian Coordinates
        double x = sphere_radius * cos(theta) * sin(phi);
        double y = sphere_radius * sin(theta) * sin(phi);
        double z = 0;

        //Since Acceleration  = GM/r^2 = M/r^2 in this system
        //Then sqrt(GM/r) = sqrt(M/r)
        //So speed is sqrt(1000/cube_width);
        //Can be deconstructed into random velocity components

        double velocity = sqrt(1000 / sphere_radius);


        struct Body new_body = {x, y, z,
                                velocity, 0, 0,
                                0, 0, 0,
                                1,
                                true};
        bodies[i] = new_body;

    }
}

void setup_random_sphere_outward(Body *bodies, int number_bodies, double sphere_radius) {
    struct Body new_body = {0, 0, 0,
                            0, 0, 0,
                            0, 0, 0,
                            1000000,
                            true};
    bodies[0] = new_body;
    double theta;
    double phi;
    for (int i = 1; i < number_bodies; i++) {

        //Random Angles
        theta = (double) 2 * M_PI * rand() / (RAND_MAX);
        phi = (double) 2 * M_PI * rand() / (RAND_MAX);

        //Polar -> Cartesian Coordinates
        double x = sphere_radius * cos(theta) * sin(phi);
        double y = sphere_radius * sin(theta) * sin(phi);
        double z = 0;

        //Since Acceleration  = GM/r^2 = M/r^2 in this system
        //Then sqrt(GM/r) = sqrt(M/r)
        //So speed is sqrt(1000/cube_width);
        //Can be deconstructed into random velocity components

        struct Body new_body = {x, y, z,
                                x, y, z, // Velocity is position (Radially outward).
                                0, 0, 0,
                                1,
                                true};
        bodies[i] = new_body;

    }
}

void setup_random_disc(Body *bodies, int number_bodies, double sphere_radius) {

    // Still picking up C++
    // Sorry about the horrible code duplication

    struct Body new_body = {0, 0, 0,
                            0, 0, 0,
                            0, 0, 0,
                            1000000,
                            true};
    bodies[0] = new_body;
    double theta;
    for (int i = 1; i < number_bodies; i++) {

        //Random Angles
        theta = (double) 2 * M_PI * rand() / (RAND_MAX);

        //Polar -> Cartesian Coordinates
        double x = sphere_radius * cos(theta);
        double y = sphere_radius * sin(theta);
        double z = 0;

        //Since Acceleration  = GM/r^2 = M/r^2 in this system
        //Then sqrt(GM/r) = sqrt(M/r)
        //So speed is sqrt(1000/cube_width);
        //Can be deconstructed into random velocity components

        double velocity = sqrt(1000 / sphere_radius);


        struct Body new_body = {x, y, z,
                                velocity, 0, 0,
                                0, 0, 0,
                                1,
                                true};
        bodies[i] = new_body;

    }
}




void save_csv_time_step(Body *bodies, int number_bodies, int counter, double loop_time, string output_name) { // Saves to output/result-counter.csv
    stringstream filename;
    filename << "output/" << output_name << "/result-" << counter << ".csv";
    ofstream out(filename.str().c_str());

    out << "x, y, z, vx, vy, vz, ax, ay, az, mass" << std::endl;
    for (int i = 0; i < number_bodies; i++) {
        if (bodies[i].active) {
            out << bodies[i].x << "," << bodies[i].y << "," << bodies[i].z << ",";
            out << bodies[i].vx << "," << bodies[i].vy << "," << bodies[i].vz << ",";
            out << bodies[i].ax << "," << bodies[i].ay << "," << bodies[i].az << ",";
            out << bodies[i].mass;
            out << endl;
            //Additional Output is for Debugging purposes in Paraview
        }

    }
    ofstream loop_file("loop_times.txt");
    loop_file << "Iteration Number:" << counter << ", Loop Time: " << loop_time << endl;

}

void reset_body_acceleration(Body *bodies, int number_bodies) {
    for (int i = 0; i < number_bodies; i++) {
        bodies[i].ax = 0.0;
        bodies[i].ay = 0.0;
        bodies[i].az = 0.0;
    }
}

void merge_bodies(Body *bodies, int i, int j) {
    cout << "Merging " << i << ", " << j;

    double momentum_x = bodies[i].vx * bodies[i].mass +
                        bodies[j].vx * bodies[j].mass;
    double momentum_y = bodies[i].vy * bodies[i].mass +
                        bodies[j].vy * bodies[j].mass;
    double momentum_z = bodies[i].vz * bodies[i].mass +
                        bodies[j].vz * bodies[j].mass;
    double total_mass = bodies[i].mass + bodies[j].mass;


    bodies[i].vx = momentum_x / total_mass;
    bodies[i].vy = momentum_y / total_mass;
    bodies[i].vz = momentum_z / total_mass;
    bodies[i].mass = total_mass;

    //Acceleration is unchanged:
    //Since bodies are close together, the rest of the
    //will have approximated them as a single body before merging
    //so only a small error is made
    bodies[j].active = false;


}

double calculate_distance(Body *bodies, int i, int j) {

    if (!bodies[i].active || !bodies[j].active) {
        cout << "Error Here!"; // Only active bodies should every reach this function
    }

    double distance = sqrt(
            (bodies[i].x - bodies[j].x) * (bodies[i].x - bodies[j].x) +
            (bodies[i].y - bodies[j].y) * (bodies[i].y - bodies[j].y) +
            (bodies[i].z - bodies[j].z) * (bodies[i].z - bodies[j].z)
    );
    return distance;
}




//Stepping Functions
void acceleration_step(Body *bodies, double MERGE_DISTANCE, int i, int j) {
    double distance = calculate_distance(bodies, i, j);
    if (distance < MERGE_DISTANCE) {
        merge_bodies(bodies, i, j);
    }
    bodies[i].ax += (bodies[j].x - bodies[i].x) * bodies[j].mass / distance / distance /
                    distance; // X Component of Force on body i exerted by body j
    bodies[i].ay += (bodies[j].y - bodies[i].y) * bodies[j].mass / distance / distance /
                    distance; // Y Component of Force on body i exerted by body j
    bodies[i].az += (bodies[j].z - bodies[i].z) * bodies[j].mass / distance / distance /
                    distance; // Z Component of Force on body i exerted by body j
}

void velocity_steps(Body *bodies, int number_bodies, double time_step_size) {
    for (int i = 0; i < number_bodies; i++) {
        if (bodies[i].active) {
            bodies[i].vx += time_step_size * bodies[i].ax;
            bodies[i].vy += time_step_size * bodies[i].ay;
            bodies[i].vz += time_step_size * bodies[i].az;
        }
    }
}

void position_steps(Body *bodies, int number_bodies, double time_step_size) {
    for (int i = 0; i < number_bodies; i++) {
        if (bodies[i].active) {
            bodies[i].x += time_step_size * bodies[i].vx;
            bodies[i].y += time_step_size * bodies[i].vy;
            bodies[i].z += time_step_size * bodies[i].vz;
        }
    }
    //cout << bodies[10].x << "," << bodies[10].y << "," << bodies[10].z << std::endl ;
}

void update_bodies(Body *bodies, int number_bodies, double MERGE_DISTANCE, double time_step_size) {
    reset_body_acceleration(bodies, number_bodies);
    for (int i = 0; i < number_bodies; i++) {
        if (bodies[i].active) {
            for (int j = 0; j < number_bodies; j++) {
                if (bodies[i].active && bodies[j].active && i != j) {
                    acceleration_step(bodies, MERGE_DISTANCE, i, j); // Acceleration exerted on i by j
                }
            }
        }
    }
    position_steps(bodies, number_bodies, time_step_size);
    velocity_steps(bodies, number_bodies, time_step_size);
}


int main(int argc, char *argv[]) {
    //Argv = [Name, Setup_type, debug, number_bodies, time_steps, time_step_size,
    // merge_distance, cube_width, sphere_radius, plot_steps, output_name]

    int setup_type = atoi(argv[1]);
    int debug = atoi(argv[2]);
    const int number_bodies = atoi(argv[3]);
    const int time_steps = atoi(argv[4]);
    const double time_step_size = atof(argv[5]);
    const double MERGE_DISTANCE = atof(argv[6]);
    const int cube_width = atoi(argv[7]);
    const int sphere_radius = atoi(argv[8]);
    const int plot_steps = atoi(argv[9]);
    const string output_name = argv[10];



    struct Body bodies[number_bodies];

    clock_t start;
    clock_t end;



    double loop_start_time = (double) clock();

    double loop_time;



    cout << "Setup Type:" << argv[1] << "\nDebug: " << argv[2] << "\n";

    if (setup_type == 1) {
        setup_random_cube(bodies, number_bodies, cube_width);
        cout << "Spherical Random Distribution \n";
    }
    else if (setup_type == 2) {
        setup_random_sphere(bodies, number_bodies, sphere_radius);
        cout << "Random Orbit\n";
    }
    else if (setup_type == 3) {
        setup_random_sphere_outward(bodies, number_bodies, sphere_radius);
        cout << "Random Orbit\n";
    }
    else if (setup_type == 4) {
        setup_random_disc(bodies, number_bodies, sphere_radius);
        cout << "Disc orbit\n";
    }
    else {
        setup_along_line(bodies, number_bodies);
        cout << "Uniform Line Distribution \n";
    }

    if (debug) {
        cout << "Debugging \n";
    }
    start = clock();
    for (int i = 0; i < time_steps; i++) {
        if (i % 100 == 0) {
            cout << "Time step " << i << endl;
        }
        update_bodies(bodies, number_bodies, MERGE_DISTANCE, time_step_size);
        if (debug) {

            double loop_end_time = (double) clock();

            loop_time = (loop_end_time - loop_start_time) / (plot_steps * CLOCKS_PER_SEC);
            if (i % plot_steps == 0) {
                save_csv_time_step(bodies, number_bodies, i, loop_time, output_name);
            }
            loop_start_time = loop_end_time;
        }
    }
    end = clock();

    double running_time = ((double) end - (double) start) / (CLOCKS_PER_SEC);
    cout << "Running Time: " << running_time << "s";
    return 0;
}
