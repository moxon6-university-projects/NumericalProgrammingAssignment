#include <fstream>
#include <sstream>
#include <math.h>
#include <iostream>
using namespace std;

struct Body{

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


const double time_step_size = 0.01;
const int time_steps = 2000000;
const int MERGE_DISTANCE = 1;
const int plot_steps = 100;


const int number_bodies = 1000;
struct Body bodies[number_bodies];

void setup_uniform(){
    for (int i = 0; i < number_bodies; i++){
        struct Body new_body = {1000 * i, 0.0, 0.0,  //Position
              0.0, 0.0, 0.0,  //Velocity
              0.0, 0.0, 0.0,  //Acceleration
              1.0,            //Mass
              true};          //Initially Active
        bodies[i] = new_body;
    }
}

/*
void setup_circular() {

}
*/

void save_csv_time_step(int counter) { // Saves to output/result-counter.csv
    stringstream filename;
    filename << "output/result-" << counter <<  ".csv";
    ofstream out( filename.str().c_str() );
    out << "x, y, z, vx, vy, vz, ax, ay, az, mass" << std::endl;
    for (int i = 0; i < number_bodies; i++) {
        if (bodies[i].active){
            out << bodies[i].x << "," << bodies[i].y << "," << bodies[i].z << ",";
            out << bodies[i].vx << "," << bodies[i].vy << "," << bodies[i].vz << ",";
            out << bodies[i].ax << "," << bodies[i].ay << "," << bodies[i].az << ",";
            out << bodies[i].mass;
            out << endl;
            //Additional Output is for Debugging purposes in Paraview
        }
    }
}

void reset_body_acceleration() {
    for (int i = 0; i < number_bodies; i++)
    {
        bodies[i].ax = 0.0;
        bodies[i].ay = 0.0;
        bodies[i].az = 0.0;
    }
}

void merge_bodies(int i, int j) {
    cout << "Merging " << i << ", " << j;

    if (bodies[i].mass < bodies[j].mass){ // Want to deactivate the body with lower mass
        int x = i;
        i = j;
        j = x;
    }

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

double calculate_distance(int i, int j) {

    if (!bodies[i].active || !bodies[j].active)
        {
            cout << "Error Here!"; // Only active bodies should every reach this function
        }

    double distance = sqrt(
            (bodies[i].x - bodies[j].x) * (bodies[i].x - bodies[j].x) +
            (bodies[i].y - bodies[j].y) * (bodies[i].y - bodies[j].y) +
            (bodies[i].z - bodies[j].z) * (bodies[i].z - bodies[j].z)
            );
    return distance;
}

void acceleration_step(int i, int j) {
    double distance = calculate_distance(i,j);
    if (distance < MERGE_DISTANCE){
        merge_bodies(i, j);
    }
    bodies[i].ax += (bodies[j].x - bodies[i].x) * bodies[j].mass / distance / distance / distance ; // X Component of Force on body i exerted by body j
    bodies[i].ay += (bodies[j].y - bodies[i].y) * bodies[j].mass / distance / distance / distance ; // Y Component of Force on body i exerted by body j
    bodies[i].az += (bodies[j].z - bodies[i].z) * bodies[j].mass / distance / distance / distance ; // Z Component of Force on body i exerted by body j
}

void velocity_steps(){
    for (int i=0; i<number_bodies; i++) {
            if (bodies[i].active){
                bodies[i].vx += time_step_size * bodies[i].ax;
                bodies[i].vy += time_step_size * bodies[i].ay;
                bodies[i].vz += time_step_size * bodies[i].az;
            }
        }
}

void position_steps() {
    for (int i=0; i<number_bodies; i++) {
        if (bodies[i].active){
            bodies[i].x += time_step_size * bodies[i].vx;
            bodies[i].y += time_step_size * bodies[i].vy;
            bodies[i].z += time_step_size * bodies[i].vz;
        }
    }
    //cout << bodies[10].x << "," << bodies[10].y << "," << bodies[10].z << std::endl ;
}

void update_bodies() {
    reset_body_acceleration();
    for (int i = 0; i < number_bodies; i++) {
        if (bodies[i].active){
            for (int j = 0; j < number_bodies; j++) {
                if (bodies[j].active && i!=j)
                {
                acceleration_step(i, j); // Acceleration exerted on i by j
                }
            }
        }
    }
    position_steps();
    velocity_steps();
}

int main(int argc, char *argv[]){

    cout <<"Setup Type:" << argv[1] << "\nDebug: " << argv[2] << "\n";


    if (argv[1] == "0"){setup_uniform();}
    else {setup_uniform();}

    bool debug;
    if (argv[2] == "0"){debug = false;}
    else{
        debug = true;
        cout << "Debugging \n";
    }




    for (int i=0; i < time_steps; i++){
        update_bodies();
        if (debug){
            if (i % plot_steps == 0)
            {
                save_csv_time_step(i);
            }
        }
    }
    return 0;
}
