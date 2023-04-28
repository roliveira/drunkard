class ParticleSystem : public RNG 
{

private:
    
    Particle part0;

    // Files

    ofstream fsummary;

    // Particles

    int dimension;
    int num_particles;

    // Grid index

    std::vector<int> grid_subindex;

    // Location

    std::vector<double> xl;
    std::vector<double> x0;
    std::vector<double> xg;

    std::vector<double> xg_mean;
    std::vector<double> xg_var;
    std::vector<double> xg_dispersion;
    std::vector<double> xg_var_old;
    std::vector<double> xg_bbox;
    std::vector<double> spacing;

    // Time

    int                 step;
    double              dt;
    double              t0;
    double              t_elapsed;
    std::vector<double> t;
    std::vector<double> tt;

    // Transitions

    std::vector<int> trans;
    int              trans_mean;
    int              trans_min;
    int              trans_max;

    // State

    std::vector<int> state;

public:

    // Constructors

    ParticleSystem(void);
    ParticleSystem(int _seed);

    // Getters
    int GetDimension(void);
    int GetNumberOfParticles(void);

    double GetPhysicalTime(int particle_index);
    double GetTransitTime (int particle_index);

    std::vector<int   > GetGridSubIndex    (int particle_index);
    std::vector<double> GetInitialPosition (int particle_index);
    std::vector<double> GetPosition        (int particle_index);
    std::vector<double> GetPeriodicPosition(int particle_index);

    std::vector<double> GetPosition(void);

    std::vector<double> GetPositionMean    (void);
    std::vector<double> GetPositionVariance(void);
    std::vector<double> GetDispersion      (void);
    std::vector<double> GetBoundingBox     (void);
    std::vector<double> GetBoundingBox0    (void);
    std::vector<double> GetBoundingBox1    (void);

    double GetDeltaTime  (void);
    double GetElapsedTime(void);
    int    GetCurrentStep(void);

    int GetTransitionsMean(void);
    int GetTransitionsMin(void);
    int GetTransitionsMax(void);

    std::vector<int> GetState(void);
    int GetState(int _index);

    // Setters

    void SetPhysicalTime         (int particle_index, double _t );
    void SetTransitTime          (int particle_index, double _tt);
    void SetPhysicalTimeIncrement(int particle_index, double _dt, double _tt);

    void SetGridSubIndex(int particle_index, int grid_subindex_x, int grid_subindex_y, int grid_subindex_z);
    void SetGridSubIndex(int particle_index, std::vector<int> grid_subindex);

    void SetPositionLocal(int particle_index, double _x, double _y, double _z);
    void SetPositionLocal(int particle_index, std::vector<double> _x);
    
    void SetPositionOrigin(int particle_index, double _x, double _y, double _z);
    void SetPositionOrigin(int particle_index, std::vector<double> _x);

    void SetPositionGlobal(int particle_index, double _x, double _y, double _z);
    void SetPositionGlobal(int particle_index, std::vector<double> _x);

    void SetDeltaTime  (double _time_delta);
    void SetInitialTime(double _time_initial);
    void SetInitialStep(int _step_initial);

    void SetInitialVariance(void);

    void SetState(std::vector<int> _states);
    void SetState(int _index, int _state);

    // Methods

    void CreateParticlePointSource(int _num_particles, std::vector<int> _shape, std::vector<double> _spacing, std::vector<double> _origin);
    void CreateParticleLineSource (int _num_particles, std::vector<int> _shape, std::vector<double> _spacing, std::vector<double> _origin);
    void CreateParticlePlaneSource(int _num_particles, std::vector<int> _shape, std::vector<double> _spacing, std::vector<double> _origin);
    void CreateParticleCubeSource (int _num_particles, std::vector<int> _shape, std::vector<double> _spacing, std::vector<double> _origin);

    void UpdateLocation(int particle_index, std::vector<double> _x, std::vector<int> _grid_subindex);

    void IncrementTime(double _t);
    void IncrementStep(void);

    void                CalculateParticleSummary(void);
    std::vector<double> CalculateMean           (std::vector<double> data);
    std::vector<double> CalculateVariance       (std::vector<double> data, std::vector<double> mean);
    std::vector<double> CalculateDispersion     (std::vector<double> variance, std::vector<double> old_variance, double delta_time);
    std::vector<double> CalculateBoundingBox    (std::vector<double> data);
    
    // I/O

    void CreateParticleSummaryFile(const char *fname);
    void WriteParticleSummary     (const char *fname);
    void WriteParticleProfile     (const char *fname);
    void WriteParticleSystem      (const char *fname);
    void WriteParticleXi          (const char *fname, double v);

};


//
// Constructors
//

inline ParticleSystem::ParticleSystem(void)
: RNG()
{}

inline ParticleSystem::ParticleSystem(int _seed) 
: RNG(_seed)
{}

//
// Getters
//

inline int ParticleSystem::GetDimension(void) 
{
    return dimension;
}

inline int ParticleSystem::GetNumberOfParticles(void)
{
    return num_particles;
}

// Time

inline double ParticleSystem::GetPhysicalTime(int particle_index) 
{
    return t.at(particle_index);
}

inline double ParticleSystem::GetTransitTime(int particle_index) 
{
    return tt.at(particle_index);
}

inline double ParticleSystem::GetElapsedTime(void) 
{
    return t_elapsed;
}

inline int ParticleSystem::GetCurrentStep(void) 
{
    return step;
}

int ParticleSystem::GetTransitionsMean(void)
{
    return trans_mean;
}

int ParticleSystem::GetTransitionsMin(void)
{
    return trans_min;
}

int ParticleSystem::GetTransitionsMax(void)
{
    return trans_max;
}

// State

inline std::vector<int> ParticleSystem::GetState(void)
{
    return state;
}

inline int ParticleSystem::GetState(int _index)
{
    return state.at(_index);
}

// Grid sub-index 

inline std::vector<int> ParticleSystem::GetGridSubIndex(int particle_index) 
{
    std::vector<int> out(dimension);
    
    for (int i = 0; i < dimension; ++i)
        out[i] = grid_subindex.at(particle_index + i * num_particles);
        
    return out;
}

// Position

inline std::vector<double> ParticleSystem::GetInitialPosition(int particle_index) 
{
    std::vector<double> out(dimension);

    for (int i = 0; i < dimension; ++i) 
        out[i] = x0.at(particle_index + i * num_particles);
    
    return out;
}

inline std::vector<double> ParticleSystem::GetPosition(int particle_index) 
{
    std::vector<double> out(dimension);

    for (int i = 0; i < dimension; ++i) 
        out[i] = xg.at(particle_index + i * num_particles);
    
    return out;
}

inline std::vector<double> ParticleSystem::GetPeriodicPosition(int particle_index) 
{
    std::vector<double> out(dimension);

    for (int i = 0; i < dimension; ++i) 
        out[i] = xl.at(particle_index + i * num_particles);
    
    return out;
}

inline std::vector<double> ParticleSystem::GetPosition(void) 
{
    return xg;
}

inline std::vector<double> ParticleSystem::GetPositionMean(void) 
{
    return xg_mean;
}

inline std::vector<double> ParticleSystem::GetPositionVariance(void) 
{
    return xg_var;
}

inline std::vector<double> ParticleSystem::GetDispersion(void) 
{
    return xg_dispersion;
}

inline std::vector<double> ParticleSystem::GetBoundingBox(void) 
{
    return xg_bbox;
}

inline std::vector<double> ParticleSystem::GetBoundingBox0(void) 
{
    return std::vector<double>(xg_bbox.begin(), xg_bbox.begin() + dimension);
}

inline std::vector<double> ParticleSystem::GetBoundingBox1(void) 
{
    return std::vector<double>(xg_bbox.begin() + dimension, xg_bbox.end());
}

//
// Setters
//

// Time

inline void ParticleSystem::SetPhysicalTime(int particle_index, double _t) 
{
    t [particle_index] = _t;
}

inline void ParticleSystem::SetTransitTime(int particle_index, double _tt) 
{
    tt[particle_index] = _tt;
}

inline void ParticleSystem::SetPhysicalTimeIncrement(int particle_index, double _dt, double _tt) 
{
    t [particle_index] += _dt;
    tt[particle_index] += _tt;
}

// Grid sub-index

inline void ParticleSystem::SetGridSubIndex(int particle_index, int _grid_subindex_x, int _grid_subindex_y, int _grid_subindex_z) 
{
    SetGridSubIndex(particle_index, {_grid_subindex_x, _grid_subindex_y, _grid_subindex_z});
}

inline void ParticleSystem::SetGridSubIndex(int particle_index, std::vector<int> _grid_subindex) 
{
    for (int i=0; i < dimension; ++i)
        grid_subindex[particle_index + i * num_particles] = _grid_subindex[i];

}

// Position

inline void ParticleSystem::SetPositionOrigin(int particle_index, double _x, double _y, double _z) 
{
    SetPositionOrigin(particle_index, {_x, _y, _z});
}

inline void ParticleSystem::SetPositionOrigin(int particle_index, std::vector<double> _x) 
{
    for (int i=0; i < dimension; ++i)
        x0[particle_index + i * num_particles] = _x[i];
}

inline void ParticleSystem::SetPositionLocal(int particle_index, double _x, double _y, double _z) 
{
    SetPositionLocal(particle_index, {_x, _y, _z});
}

inline void ParticleSystem::SetPositionLocal(int particle_index, std::vector<double> _x) 
{
    for (int i=0; i < dimension; ++i)
        xl[particle_index + i * num_particles] = _x[i];
}

inline void ParticleSystem::SetPositionGlobal(int particle_index, double _x, double _y, double _z) 
{
    SetPositionGlobal(particle_index, {_x, _y, _z});
}

inline void ParticleSystem::SetPositionGlobal(int particle_index, std::vector<double> _x) 
{
    for (int i=0; i < dimension; ++i)
        xg[particle_index + i * num_particles] = _x[i];
}

inline void ParticleSystem::SetDeltaTime(double _time_delta) 
{
    dt = _time_delta;
}

inline void ParticleSystem::SetInitialTime(double _time_initial)
{
    t0        = _time_initial;
    t_elapsed = _time_initial;
}

inline void ParticleSystem::SetInitialStep(int _step_initial)
{
    step = _step_initial;
}

inline void ParticleSystem::SetInitialVariance(void)
{
    xg_var_old = std::vector<double>(dimension, 0.0);
}

// State

inline void ParticleSystem::SetState(std::vector<int> _states)
{
    state = _states;
}

inline void ParticleSystem::SetState(int _index, int _state)
{
    state[_index] = _state;
}

//
// Methods
//

// Create particles

inline void ParticleSystem::CreateParticlePointSource(int _num_particles, std::vector<int> _shape, std::vector<double> _spacing, std::vector<double> _origin) 
{
    dimension     = _shape.size();
    spacing       = _spacing;
    num_particles = _num_particles;
   
    state = std::vector<int>(num_particles, 1);

    trans = std::vector<int>(num_particles, 0);
    grid_subindex = std::vector<int>(dimension * num_particles);

    xl = std::vector<double>(dimension * num_particles);
    x0 = std::vector<double>(dimension * num_particles);
    xg = std::vector<double>(dimension * num_particles);

    t  = std::vector<double>(num_particles);
    tt = std::vector<double>(num_particles);

    for(int ip = 0; ip < num_particles; ++ip) 
    {
        std::vector<int   > temp_i(dimension);
        std::vector<double> temp_x(dimension);

        for (int i = 0; i < dimension; ++i) 
        {
            temp_i[i] = static_cast<int>(_shape[i] / 2.0);
            temp_x[i] = _origin[i] + static_cast<double>(temp_i[i]) * _spacing[i];
        }

        SetGridSubIndex  (ip, temp_i);
        SetPositionOrigin(ip, temp_x);
        SetPositionLocal (ip, temp_x);
        SetPositionGlobal(ip, temp_x);

        SetPhysicalTime(ip, 0.0);
        SetTransitTime (ip, 0.0);
    }    
}

inline void ParticleSystem::CreateParticleLineSource(int _num_particles, std::vector<int> _shape, std::vector<double> _spacing, std::vector<double> _origin) 
{
    dimension     = _shape.size();
    spacing       = _spacing;
    num_particles = _num_particles;
    
    state = std::vector<int>(num_particles, 1);

    std::vector<std::uniform_int_distribution<int>> dist(dimension);
    for (int i=0; i < dimension; ++i)
        dist[i] = std::uniform_int_distribution<int>(1, _shape[i] - 2);

    trans = std::vector<int>(num_particles, 0);
    grid_subindex = std::vector<int>(dimension * num_particles);

    xl = std::vector<double>(dimension * num_particles);
    x0 = std::vector<double>(dimension * num_particles);
    xg = std::vector<double>(dimension * num_particles);

    t  = std::vector<double>(num_particles);
    tt = std::vector<double>(num_particles);

    for(int ip = 0; ip < num_particles; ++ip) {
        std::vector<int   > temp_i(dimension);
        std::vector<double> temp_x(dimension);

        for (int i = 0; i < dimension; ++i) {
            if (i == 1) temp_i[i] = dist[i](this->_generator);
            else        temp_i[i] = _shape[i] / 2;

            temp_x[i] = _origin[i] + static_cast<double>(temp_i[i]) * _spacing[i];
        }

        SetGridSubIndex  (ip, temp_i);
        SetPositionOrigin(ip, temp_x);
        SetPositionLocal (ip, temp_x);
        SetPositionGlobal(ip, temp_x);

        SetPhysicalTime(ip, 0.0);
        SetTransitTime (ip, 0.0);
    }
}

inline void ParticleSystem::CreateParticlePlaneSource(int _num_particles, std::vector<int> _shape, std::vector<double> _spacing, std::vector<double> _origin) 
{
    dimension     = _shape.size();
    spacing       = _spacing;
    num_particles = _num_particles;
    
    state = std::vector<int>(num_particles, 1);

    std::vector<std::uniform_int_distribution<int>> dist(dimension);
    for (int i=0; i < dimension; ++i)
        dist[i] = std::uniform_int_distribution<int>(1, _shape[i] - 2);

    trans = std::vector<int>(num_particles, 0);
    grid_subindex = std::vector<int>(dimension * num_particles);

    xl = std::vector<double>(dimension * num_particles);
    x0 = std::vector<double>(dimension * num_particles);
    xg = std::vector<double>(dimension * num_particles);

    t  = std::vector<double>(num_particles);
    tt = std::vector<double>(num_particles);

    for(int ip = 0; ip < num_particles; ++ip) {
        std::vector<int   > temp_i(dimension);
        std::vector<double> temp_x(dimension);

        for (int i = 0; i < dimension; ++i) {
            if (i == 0) temp_i[i] = 1;
            else        temp_i[i] = dist[i](this->_generator);

            temp_x[i] = _origin[i] + static_cast<double>(temp_i[i]) * _spacing[i];
        }

        SetGridSubIndex  (ip, temp_i);
        SetPositionOrigin(ip, temp_x);
        SetPositionLocal (ip, temp_x);
        SetPositionGlobal(ip, temp_x);

        SetPhysicalTime(ip, 0.0);
        SetTransitTime (ip, 0.0);
    }
}

inline void ParticleSystem::CreateParticleCubeSource(int _num_particles, std::vector<int> _shape, std::vector<double> _spacing, std::vector<double> _origin) 
{
    dimension     = _shape.size();
    spacing       = _spacing;
    num_particles = _num_particles;
    
    state = std::vector<int>(num_particles, 1);

    std::vector<std::uniform_int_distribution<int>> dist(dimension);
    for (int i=0; i < dimension; ++i)
        dist[i] = std::uniform_int_distribution<int>(1, _shape[i] - 2);

    trans = std::vector<int>(num_particles, 0);
    grid_subindex = std::vector<int>(dimension * num_particles);

    xl = std::vector<double>(dimension * num_particles);
    x0 = std::vector<double>(dimension * num_particles);
    xg = std::vector<double>(dimension * num_particles);

    t  = std::vector<double>(num_particles);
    tt = std::vector<double>(num_particles);

    for(int ip = 0; ip < num_particles; ++ip) {
        std::vector<int   > temp_i(dimension);
        std::vector<double> temp_x(dimension);

        for (int i = 0; i < dimension; ++i) {
            temp_i[i] = dist[i](this->_generator);
            temp_x[i] = _origin[i] + static_cast<double>(temp_i[i]) * _spacing[i];
        }

        SetGridSubIndex  (ip, temp_i);
        SetPositionOrigin(ip, temp_x);
        SetPositionLocal (ip, temp_x);
        SetPositionGlobal(ip, temp_x);

        SetPhysicalTime(ip, 0.0);
        SetTransitTime (ip, 0.0);
    }
}

inline void ParticleSystem::UpdateLocation(int particle_index, std::vector<double> _x, std::vector<int> _grid_subindex) 
{
    SetGridSubIndex  (particle_index, _grid_subindex);
    SetPositionLocal (particle_index, _x            );
    SetPositionGlobal(particle_index, _x            );
}

inline void ParticleSystem::IncrementTime(double _t) 
{
    t_elapsed += _t;
}

inline void ParticleSystem::IncrementStep(void) 
{
    step += 1;
}

inline void ParticleSystem::CalculateParticleSummary(void) 
{
    xg_mean       = CalculateMean(xg);
    xg_var        = CalculateVariance(xg, xg_mean);
    xg_dispersion = CalculateDispersion(xg_var, xg_var_old, dt);
    xg_bbox       = CalculateBoundingBox(xg);
    
    trans_mean    = accumulate(trans.begin(), trans.end(), 0.0) / static_cast<float>(num_particles);

    std::pair<std::vector<int>::iterator, std::vector<int>::iterator> val = 
            std::minmax_element(trans.begin(), trans.end());
    trans_min = *val.first;
    trans_max = *val.second;

    xg_var_old = xg_var;
}

inline std::vector<double> ParticleSystem::CalculateMean(std::vector<double> data) 
{
    std::vector<double> mean(dimension, 0.0);

    for (int i = 0; i < dimension; ++i) {
        for (int ip = 0; ip < num_particles; ++ip)
            mean[i] += data[ip + i * num_particles];
            
        mean[i] /= static_cast<double>(num_particles);
    }

    return mean;
}

inline std::vector<double> ParticleSystem::CalculateVariance(std::vector<double> data, std::vector<double> mean) 
{
    std::vector<double> var(dimension);

    for (int i = 0; i < dimension; ++i) {
        double temp_var = 0;
        
        for (int ip = 0; ip < num_particles; ++ip) {
            temp_var += std::pow(data[ip + i * num_particles] - mean[i], 2);
        }
        
        var[i] = temp_var / static_cast<double>(num_particles);
    }

    return var;
}

inline std::vector<double> ParticleSystem::CalculateDispersion(std::vector<double> variance, std::vector<double> old_variance, double delta_time) 
{
    std::vector<double> disp(dimension);

    for (int i = 0; i < dimension; ++i) {
        disp[i] = (variance[i] - old_variance[i]) / delta_time;
    }

    return disp;
}

inline std::vector<double> ParticleSystem::CalculateBoundingBox(std::vector<double> data) 
{
    std::vector<double> bbox(2 * dimension);

    for(int i = 0; i < dimension; ++i) 
    {
        std::pair<std::vector<double>::iterator, std::vector<double>::iterator> val = 
            std::minmax_element(data.begin() + i * num_particles, data.begin() + (i + 1) * num_particles);
            
        bbox[i       ] = *val.first;
        bbox[i + dimension] = *val.second;
    }

    return bbox;
}

//
// I/O
//

// Writers

inline void ParticleSystem::WriteParticleSystem(const char *fname) {
    
    // Point data collection

    vtkSmartPointer<vtkDataArrayCollection> point_collection =
        vtkSmartPointer<vtkDataArrayCollection>::New();

    // Field data collection

    vtkSmartPointer<vtkDataArrayCollection> field_collection = 
        vtkSmartPointer<vtkDataArrayCollection>::New();

    // Point data arrays

    vtkSmartPointer<vtkFloatArray> xl_array =
        vtkSmartPointer<vtkFloatArray>::New();
    xl_array->SetName("xl");
    xl_array->SetNumberOfComponents(dimension);
    xl_array->Allocate(num_particles);
    point_collection->AddItem(xl_array);

    vtkSmartPointer<vtkFloatArray> x0_array =
        vtkSmartPointer<vtkFloatArray>::New();
    x0_array->SetName("x0");
    x0_array->SetNumberOfComponents(dimension);
    x0_array->Allocate(num_particles);
    point_collection->AddItem(x0_array);

    vtkSmartPointer<vtkFloatArray> xg_array =
        vtkSmartPointer<vtkFloatArray>::New();
    xg_array->SetName("xg");
    xg_array->SetNumberOfComponents(dimension);
    xg_array->Allocate(num_particles);
    point_collection->AddItem(xg_array);

    vtkSmartPointer<vtkFloatArray> t = 
        vtkSmartPointer<vtkFloatArray>::New();
    t->SetName("t");
    t->SetNumberOfComponents(1);
    t->Allocate(num_particles);
    point_collection->AddItem(t);

    vtkSmartPointer<vtkFloatArray> tt = 
        vtkSmartPointer<vtkFloatArray>::New();
    tt->SetName("tt");
    tt->SetNumberOfComponents(1);
    tt->Allocate(num_particles);
    point_collection->AddItem(tt);

    vtkSmartPointer<vtkIntArray> state = 
        vtkSmartPointer<vtkIntArray>::New();
    state->SetName("state");
    state->SetNumberOfComponents(1);
    state->Allocate(num_particles);
    point_collection->AddItem(state);

    // Field data arrays

    vtkSmartPointer<vtkIntArray> step_array = 
        vtkSmartPointer<vtkIntArray>::New();
    step_array->SetName("step");
    step_array->SetNumberOfComponents(1);
    field_collection->AddItem(step_array);

    vtkSmartPointer<vtkFloatArray> t0_array = 
        vtkSmartPointer<vtkFloatArray>::New();
    t0_array->SetName("t0");
    t0_array->SetNumberOfComponents(1);
    field_collection->AddItem(t0_array);

    vtkSmartPointer<vtkFloatArray> t_array = 
        vtkSmartPointer<vtkFloatArray>::New();
    t_array->SetName("t");
    t_array->SetNumberOfComponents(1);
    field_collection->AddItem(t_array);

    vtkSmartPointer<vtkFloatArray> xg_mean_array = 
        vtkSmartPointer<vtkFloatArray>::New();
    xg_mean_array->SetName("xg_mean");
    xg_mean_array->SetNumberOfComponents(1);
    field_collection->AddItem(xg_mean_array);

    vtkSmartPointer<vtkFloatArray> xg_var_array = 
        vtkSmartPointer<vtkFloatArray>::New();
    xg_var_array->SetName("xg_var");
    xg_var_array->SetNumberOfComponents(1);
    field_collection->AddItem(xg_var_array);

    vtkSmartPointer<vtkFloatArray> xg_dispersion_array = 
        vtkSmartPointer<vtkFloatArray>::New();
    xg_dispersion_array->SetName("xg_dispersion");
    xg_dispersion_array->SetNumberOfComponents(1);
    field_collection->AddItem(xg_dispersion_array);

    vtkSmartPointer<vtkFloatArray> xg_bbox_array = 
        vtkSmartPointer<vtkFloatArray>::New();
    xg_bbox_array->SetName("xg_bbox");
    xg_bbox_array->SetNumberOfComponents(2 * dimension);
    field_collection->AddItem(xg_bbox_array);

    // Points

    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
    points->Allocate(num_particles);

    // Cells

    vtkSmartPointer<vtkCellArray> verts =
        vtkSmartPointer<vtkCellArray>::New();

    // Set data array values from nodes
    int count = 0;

    for (int ip = 0; ip < num_particles; ++ip) {
        vtkIdType vtkid[1] = { count };  // Confused why you can't use ints or
                                         // unsigned ints in some places?
                                         // http://www.vtk.org/Wiki/VTK/Tutorials/VtkIdType

        points->InsertNextPoint(GetPosition(ip).data());

        // points->InsertNextPoint(xg[ip], xg[ip + num_particles], xg[ip + 2 * num_particles]);

        verts->InsertNextCell(1, vtkid);
        
        xl_array->InsertNextTuple(GetPeriodicPosition(ip).data());
        x0_array->InsertNextTuple(GetInitialPosition(ip).data());
        xg_array->InsertNextTuple(GetPosition(ip).data());

        t->InsertNextValue(GetPhysicalTime(ip));
        tt->InsertNextValue(GetTransitTime(ip));

        state->InsertNextValue(GetState(ip));

        ++count;
    }

    t_array->InsertNextValue(t_elapsed);
    t0_array->InsertNextValue(t0);
    step_array->InsertNextValue(step);
    xg_mean_array->InsertNextTuple(xg_mean.data());
    xg_var_array->InsertNextTuple(xg_var.data());
    xg_dispersion_array->InsertNextTuple(xg_dispersion.data());
    xg_bbox_array->InsertNextTuple(xg_bbox.data());

    // PolyData dataset

    vtkSmartPointer<vtkPolyData> polydata =
        vtkSmartPointer<vtkPolyData>::New();

    polydata->SetPoints(points);
    polydata->SetVerts (verts);

    for (int i = 0; i < point_collection->GetNumberOfItems(); ++i)
        polydata->GetPointData()->AddArray(point_collection->GetItem(i));

    for (int i = 0; i < field_collection->GetNumberOfItems(); ++i)
        polydata->GetFieldData()->AddArray(field_collection->GetItem(i));   

    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(fname);
    writer->SetInputData(polydata);
    writer->Write();
}

inline void ParticleSystem::CreateParticleSummaryFile(const char *fname) {
    fsummary.open(fname, std::ofstream::out);

    fsummary << std::right 
             << /*setw(15) <<*/ "step" << ","
             << /*setw(15) <<*/ "time" << ",";

    for (int i=0; i<dimension;   ++i) fsummary << /*setw(15) <<*/ "x" << i << "_min"  << ",";
    for (int i=0; i<dimension;   ++i) fsummary << /*setw(15) <<*/ "x" << i << "_max"  << ",";
    for (int i=0; i<dimension;   ++i) fsummary << /*setw(15) <<*/ "x" << i << "_mean" << ",";
    for (int i=0; i<dimension;   ++i) fsummary << /*setw(15) <<*/ "x" << i << "_var"  << ",";
    for (int i=0; i<dimension-1; ++i) fsummary << /*setw(15) <<*/ "x" << i << "_dl"   << ",";

    if (dimension-1 > 0) fsummary << /*setw(15) <<*/ "x" << dimension-1 << "_dl";

    fsummary << std::endl;
    fsummary.close();
}

inline void ParticleSystem::WriteParticleSummary(const char *fname) {
    fsummary.open(fname, std::ofstream::out | std::ofstream::app);

    std::vector<double> _bb = GetBoundingBox();
    std::vector<double> _xm = GetPositionMean();
    std::vector<double> _xv = GetPositionVariance();
    std::vector<double> _dl = GetDispersion();

    fsummary << std::right 
             << /*setw(15) <<*/ step      << ","
             << /*setw(15) <<*/ t_elapsed << ",";

    for (int i=0; i<2*dimension;   ++i) fsummary << /*setw(15) <<*/ _bb[i] << ",";
    for (int i=0; i<dimension;     ++i) fsummary << /*setw(15) <<*/ _xm[i] << ",";
    for (int i=0; i<dimension;     ++i) fsummary << /*setw(15) <<*/ _xv[i] << ",";
    for (int i=0; i<dimension-1;   ++i) fsummary << /*setw(15) <<*/ _dl[i] << ",";
    
    if (dimension-1 > 0) fsummary << /*setw(15) <<*/ _dl[dimension-1] << " ";

    fsummary << std::endl;
    fsummary.close();
}

// inline void ParticleSystem::WriteParticleProfile(const char *fname) {
//     if (t_elapsed > 0.0) {
//         histogram(
//             std::vector<double>(xg.begin(), xg.begin() + num_particles),
//             xg_bbox[0], 
//             xg_bbox[dimension], 
//             static_cast<int>((xg_bbox[dimension] - xg_bbox[0]) / spacing[0]), 
//             fname
//         );
//     }
// }

// inline void ParticleSystem::WriteParticleXi(const char *fname, double v) {
//     if (t_elapsed > 0.0) {
//         std::vector<double> xi = std::vector<double>(xg.begin(), xg.begin() + num_particles);
//         for(int i = 0; i < num_particles; ++i) xi[i] = (xi[i] - x0[i]) / (v * t_elapsed);
//         histogram(xi, 0.0, 5.0, 20, fname);
//     }
// }

// Streams 

inline std::ostream &operator<<(std::ostream &os, const std::shared_ptr<ParticleSystem> &ps) 
{
    os << message("Current iteration",      ps->GetCurrentStep()      );
    os << message("Number of particles",    ps->GetNumberOfParticles());
    os << message("Elapsed time [s]",       ps->GetElapsedTime()      );
    os << message("Transitions min",        ps->GetTransitionsMin()   );
    os << message("Transitions mean",       ps->GetTransitionsMean()  );
    os << message("Transitions max",        ps->GetTransitionsMax()   );
    os << message("Position min [m]",       ps->GetBoundingBox0()     );
    os << message("Position mean [m]",      ps->GetPositionMean()     );
    os << message("Position max [m]",       ps->GetBoundingBox1()     );
    os << message("Position variance [m2]", ps->GetPositionVariance() );
    os << message("Dispersion [m]",         ps->GetDispersion()       );
    
    return os;
}