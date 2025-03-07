
__global__ void matrix_vector_kernel(
    double * lhs, 
    double * map1,
    double * map2,
    double * matrix,
    double * x
)
{
    lhs(map1(df,cell)+k) = lhs(map1(df,cell)+k) + matrix(ik,df,df2)*x(map2(df2,cell)+k);
}