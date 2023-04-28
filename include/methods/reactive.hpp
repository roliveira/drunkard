
#ifndef METHODS_REACTIVE_HPP_
#define METHODS_REACTIVE_HPP_


class Reactive {
public:

    double kr;
    double ceq;
    double dt;

    Reactive(double kr, double ceq, double dt);

    double dcdt(double kr, double ceq, double t, double c);
    double dc  (double kr, double ceq, double t, double c);
};


Reactive::Reactive(double kr, double ceq, double dt) : kr(kr), ceq(ceq), dt(dt) {
}

double Reactive::dcdt(double kr, double ceq, double t, double c) {
    return kr*(ceq-c);
}

double Reactive::dc(double kr, double ceq, double t, double c) {

    double t_temp = 0.0;
    double c_temp = c;

    while (t_temp < t) {
        double k1 = dt*this->dcdt(kr, ceq, t_temp             , c       );
        double k2 = dt*this->dcdt(kr, ceq, t_temp+this->dt/2.0, c+k1/2.0);
        double k3 = dt*this->dcdt(kr, ceq, t_temp+this->dt/2.0, c+k2/2.0);
        double k4 = dt*this->dcdt(kr, ceq, t_temp+this->dt    , c+k3    );

        c_temp += (k1+2*k2+2*k3+k4)/6.0;
        t_temp += dt;
    }

    return c_temp-c;
}


#endif  // METHODS_REACTIVE_HPP_
