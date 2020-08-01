# Project Topic
## Optimization for Polymer Predictors for Drug Delivery Systems, using Machine Learning Approaches

### MOTIVATION
Polymer designing for a particular drug (usually a poorly soluble drug) is a very tedious and time-consuming process. 
Hence, the entire challenge for my work is to reduce the amount of effort required for experimentation. I used several 
Machine Learning approaches to develop correlations between the predictors, which could make the experiments easier and 
time saving. Finally, we were able to propose a set of protocols, which could be followed for choosing the best polymer 
for a particular drug.

<img width="1438" alt="Screen Shot 2020-08-02 at 12 23 50 AM" src="https://user-images.githubusercontent.com/66521525/89108646-abee8400-d457-11ea-8da0-682e8a407b8c.png">

This repository consists of the proposed code to be followed to clean the raw dataset and subsequently shorlisting the best 
subset of polymers for drug, NIFEDIPINE.

### IMPLEMENTATION is divided into the following parts :-
#### 1) Principal Component Analysis (PCA) - 
Part of Data Mining and Data Visualization. The entire polymer data was studied in order to look at :–
* How well the data variation is being covered by the predictors.
* Which all predictors are more important than the others.
* What is the correlation between thedifferent predictors.
* Can removing some predictors from the dataset actually make a better model?

#### 2) Multivariate Linear Regression Model
Entire work used a multiple linear regression model.

#### 3) Best Subset Selection Method - 
Part of Feature Selection and Data Regularization. Chose the best subset for each predictor as the response.

### CONCLUSION
Proposed a set of protocols, which could be followed for the prediction of Tg and Sol in all the 3 scenarios :–
* If required data is available.
* If required data isn’t available, but some data is available.
* If required data is not available and some other data required to predict the relevant data is also not available.
