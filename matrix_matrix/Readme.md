# Matrix matrix kernel 

| Paradigm | Time($\mu s$) | Notes
--- | --- | ---
| cuda | 945.616 | Naive implememntation + thread ordering (x,y) 
| openmp | 972 | Naive implememntation
| openmp | 1011.4 | Naive impl. + scheduling manual (as in cuda version)
| cuda | 235.295 | Naive impl. + shared_mem  + thread ordering (x,y)
| openmp | 634.644 | Naive impl. + shared_mem
| openmp | 560.093 | Naive impl. + shared_mem + loop inside parallel region