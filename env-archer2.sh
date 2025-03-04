

if [[ -z "${OPENMP_OFFLOAD_COURSE_ROOT}" ]]; then
    
    module load cpe/23.09
    module load PrgEnv-cray
    module load rocm
    module load craype-accel-amd-gfx90a
    module load craype-x86-milan
    module load cray-python


    OPENMP_OFFLOAD_COURSE_ROOT="${HOME/home/work}/.openmp_offload_course"


    if [ ! -d "$OPENMP_OFFLOAD_COURSE_ROOT/python_env" ]; then
        echo "Python testing virtual environment not found, creating it now..."
        python3 -m venv "$OPENMP_OFFLOAD_COURSE_ROOT/python_env"
        source "$OPENMP_OFFLOAD_COURSE_ROOT/python_env/bin/activate"
        pip install -r requirements.txt
        source "$OPENMP_OFFLOAD_COURSE_ROOT/python_env/bin/activate"
    else
        source "$OPENMP_OFFLOAD_COURSE_ROOT/python_env/bin/activate"
    fi
fi