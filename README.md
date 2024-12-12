


# SphericalRotationDimensionReduction (SRCA)

**Content**
This is the code repository for the research publication "Spherical Rotation Dimension Reduction with Geometric Loss Functions (SRCA)" by [Didong Li](), [Jeremy E. Purvis](https://www.med.unc.edu/genetics/purvislab/) and [Hengrui Luo](https://hrluo.github.io/SphericalRotationDimensionReduction). 
The manuscript of this paper can be accessed at https://www.jmlr.org/papers/v25/23-0547.html or https://arxiv.org/abs/2204.10975. 

 - The MATLAB code we store in root folder contains a full pipeline for two methods:
 	- Spherical rotation (SRCA) dimension reduction.
	- Spherlets approximation dimension reduction, see [Efficient Manifold and SubspaceApproximations with Spherelets](https://arxiv.org/abs/1706.08263). 
	- Principal component analysis for comparison convenience. 
 - The folder stores the dataset we used in our paper, with example code (i.e.,`BasicExample_*` for synthetic examples; `Example_*` for real data analysis) for both basic (synthetic) example and the data example (i.e.,`Example_*`). 
 - The file `randvonMisesFisherm.m` belongs to Dr. [Sungkyu Jung](https://www.stat.pitt.edu/sungkyu/oldSoftwarePage.html), it is only used for generating synthetic data and is *NOT* used in the SRCA methodology. 
 - The file `Isomap.m` belongs to Tenenbaum, de Silva, and Langford (2000). 

**Abstract**
Modern datasets witness high-dimensionality and nontrivial geometries of spaces they live in. It would be helpful in data analysis that we can reduce the dimensionality while retaining the geometric structure of the dataset. 
Motivated by this observation, 
we want to incorporate geometrical information in the task of dimension reduction. Following this thought, we propose a general method where the low-dimensional manifold is in the shape of ellipsoids in this paper. Our Spherical Rotation Component Analysis (SRCA) is a dimension reduction method that uses spheres or ellipsoids, to approximate a low-dimensional manifold. This method not only generalizes the Spherelets method in terms of theories and algorithms and presents a comprehensive comparison of our method, as an optimization problem with theoretical guarantee and also as a structural preserving low-rank representation of data. Results relative to state-of-the-art competitors show considerable gains in ability to accurately approximate the subspace with fewer components and better structural preserving. In addition, we have pointed out that this method is a specific incarnation of a grander idea of using a geometrically induced loss function in dimension reduction tasks.

**Citation**
We provided MATLAB code for reproducible and experimental purposes under [LICENSE](https://github.com/hrluo/SphericalRotationDimensionReduction).
Please cite our paper using following BibTeX item:

	@article{luo2024spherical,
	  title={Spherical rotation dimension reduction with geometric loss functions},
	  author={Luo, Hengrui and Purvis, Jeremy E and Li, Didong},
	  journal={Journal of Machine Learning Research},
	  volume={25},
	  number={175},
	  pages={1--55},
	  year={2024}
	}

Thank you again for the interest and please reach out if you have further questions.
