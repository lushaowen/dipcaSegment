Segmentation of Multivariate Industrial Time Series Data Based on Dynamic Latent Variable Predictability
====

# dipcaSegment

This project implements an time series segmentation algorithm introduced in the following paper:

>>Shaowen Lu and Shuyu Huang, "Segmentation of multivariate industrial time series data based on dynamic latent variable predictability," IEEE Access, vol. 8, pp. 112092–112103, 2020. 
>>[open access pdf](https://ieeexplore.ieee.org/document/9116988)

The proposed algorithm recursively merges neighborhood subsequences through a heuristic bottom-up procedure. The cost of merging is defined as the mutual predictability of the subsequence models which are constructed using the [DiPCA algorithm](https://www.sciencedirect.com/science/article/pii/S095915241730094X). Then, the goal becomes finding the segmentation scheme which minimizes the averaged dynamic prediction errors of each subsequence model.

Please cite the paper if using the code:
```
@ARTICLE{Lu2020Access,
  author = {S. {Lu} and S. {Huang}},
  title = {Segmentation of Multivariate Industrial Time Series Data Based on
	Dynamic Latent Variable Predictability},
  journal = {IEEE Access},
  year = {2020},
  volume = {8},
  pages = {112092-112103},
  doi = {10.1109/ACCESS.2020.3002257},
  issn = {2169-3536}
  }
```

# usage
(1) run demo example: "demo_DiPCAseg_lasso.m"
(2) run "synthetical_TS_generate.m" to generate new data

# license
MIT © Shaowen Lu
