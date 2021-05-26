# Diffusion Map of the Dietary Niche

## From the N-D Hutchinsonian Niche to the Niche Manifold

We will now address how we will we compare the foraging behavior of organisms operating within complex seasonal environments, and relate these behaviors to measures of fitness in a framework that is investigative and potentially predictive.
Hutchinson famously proposed an N-dimensional hypervolume to describe an organism's niche, however we must suppose that most of this potential volume cannot be realized due to physical, behavioral, or environmental constraints. 
We can therefore imagine a process by which the potential space can be observed apart from that which is nonphysical or unrealizable for a particular suite of organisms to gain a clearer understanding of species' roles within a larger community.
In the context of the dietary niche, mechanistic foraging models provide a means by which to quantify the potential niche space given a range of constraints, as the consequences of a large range of foraging strategies can be systematically assessed -- in the albeit simplified manner of mathematical models.
If the constraints captured by the foraging model play a central role in shaping the foraging behaviors of consumers in natural systems, the niche-space enumerated by the model can be used as a lens by which the larger community can be understood.
However, even a simple foraging model can result in an enormous complexity of foraging strategies, and a central challenge in understanding both model and natural ecological systems resides in excising low-dimensional relationships from high-dimensional data, particularly when such data are generated from interactions with nonlinear dependencies.

Here we propose to devise a series of foraging models that range from simple to complex, to explore and enumerate the dietary consequences of a large range of foraging strategies.
We will then use this range of strategies to construct a low-dimensional realization of a high-dimensional niche space using diffusion mapping techniques.
Diffusion mapping techniques constitute a class of data analytics capable of reconstructing nonlinear relationships governing high-dimensional datasets when the generative processes underlying the data is unknown.
We will use this perspective to reconstruct a low-dimensional embedding of consumer strategies, or *strategy-niche manifold*, capturing the structural similarities of a range of foraging strategies.
Understanding how different strategies relate to one another when the mechanics are simple and well-defined provides a null expectation by which to classify and evaluate the more complex strategies observed among consumers in natural systems.
Our procedure thus consists of four parts: 1) build a class of foraging models to systematically capture a range of foraging strategies given known physical and biological constraints; 2) reconstruct a low-dimensional embedding of the range of modeled foraging behaviors to capture the associated *strategy-niche manifold*; 3) measure the diets of local consumers and assess their ecological roles by their proximity to modeled strategies on the niche manifold. 
<!-- Next, we will use this framework to 4) examine whether aspects of consumer fitness are the *niche manifold*. -->

Because both the underlying generative models and resulting foraging strategies of model systems within the embedding are known, in addition to assessment, we can use the constructed niche manifold in a predictive capacity.
For instance, the foraging constraints that give rise to different strategies also allow us to assess the fitness consequences of those strategies.
By relating the fitness consequences of simulated strategies occupying different locations along the niche manifold, we can associate empirically measured strategies with modeled values nearby.
The resultant fitness landscape assessed along the niche manifold can thus be used to evaluate and ultimately predict life history characteristics of consumers foraging in natural systems.
Thus, by incorporating foraging data and estimates of fitness from Sevilleta consumers alongside those simulated within the context of the niche manifold, we will directly assess the predictive value of the manifold niche concept in a natural community.


## A simple foraging null model to define consumer strategies

We present a minimal consumer foraging model to illustrate how our framework can be used to uncover a simulated consumer's niche manifold -- which we will identify using diffusion mapping techniques -- and treat this as a null expectation by which we will subsequently evaluate empirical consumer strategies.
We note that mechanistic foraging models of arbitrary complexity could be used to establish this null expectation, and while we will identify in what direction we aim to build additional complexity, a simpler framework best suits our illustration of the core approach.

We simulate a consumer foraging in a seasonal environment given a set of resource functional groups (see xx) with differing spatial distributions during fall (non-monsoonal) and spring (monsoonal) seasons, such that mu_si is the mean encounter rate where s denotes season and i=1...N functional groups.
A consumer of mass M forages within this landscape, targeting a particular functional group with weight tau.
The targeting weight tau indicates that for each consumer-resource interaction, the consumer will find and acquire its *targeted* resource with probability tau, and target the *closest* resource group (regardless of preference) with probably 1-tau.
Once a consumer-resource interaction is drawn, the consumer travels the distance to the resource with velocity v(M), and assimilates both bulk energy (kJ) and nitrogen content based on the resource's energy density and nitrogen concentration, respectively.
Consumer-resource interactions continue until a predetermined time threshold for the foraging bout, tmax, is reached, whereupon the consumer ceases its foraging activity for the day.
We track 300 days of consumer foraging, with 100 days committed to each season in a `fall-spring-fall` cycle.


## Exploring consumer niche-space with diffusion mapping

The simplistic foraging model described above allows us to track both temporal changes in diet as well as daily energetic gains as a function of consumer resource targeting strategies.
Where a given consumer targeting strategy is characterized by the temporal sequence of the proportional contribution of each resource group to the consumer's diet, we then use a diffusion mapping approach to compare alternative consumer targeting strategies.
In the example we present here, there are 7 resource groups and we explore targeting strategies ranging from tau=0 (a consumer that always targets the nearest resource) to tau=1 (a consumer that always targets one of the 7 resource groups regardless of distance).
Resource groups include C3/C4 perennial/annual shrubs, forbs, and grasses with seasonal encounter rates scaled to seasonal densities.
Accordingly a single consumer strategy is described by the targeting of a specific resource group with weight tau, and exists as a single point within the diffusion space. 

After averaging biweekly dietary vectors of each targeting strategy across n=500 replicate foraging simulations, and following the principles of diffusion mapping, we first establish a similarity matrix across targeting strategy pairs.
A given targeting strategy can be represented by a matrix of biweekly proportional contribution averages (with rows being resource groups, and columns being biweekly averages), and we calculate pairwise similarity between each matrix pair using Jaccard distance.
The similarity matrix can be treated as an Adjacency matrix A - defining how nodes in a network are linked - where strategies (nodes) are linked together if they have a similarity greater than a particular value, and are not linked if they fall below that value.
In this case, each targeting strategy is linked to the k=10 most similar targeting strategies, and we note that our results are not particularly sensitive to the choice of k. 
We then imagine a diffusive process taking place on this `strategy network', where the diffusive modes can be used to construct a consumer strategy space where alternative strategies can be directly compared to one another.
To obtain these diffusive modes, we transform the strategy matrix into a Laplacian matrix, such that L = A - D, where D is the diagonal matrix of A. 
The eigenvectors of L (v_0 to v_n) thus provide the modes of the diffusive process, which can be scaled to the Laplacian eigenvalues (lambda_0 to lambda_n). 
From the the n-dimensional diffusion space and following Fahimipour et al. (2020), we finally construct a 2-dimensional embedding of the consumer strategy-niche manifold, to permit easier visualization.
Strategies that share greater temporal similarity with each other will thus fall closer together within the diffusion eigenspace, and form the basis by which alternative empirical strategies can then be assessed.


## Preliminary results of the consumer strategy-niche manifold

The 2-D embedding of the consumer strategy-niche manifold depicts the array of simulated consumer targeting strategies as points, where point color denotes the resource group targeted and opacity denotes targeting weight (increased opacity means that the consumer targets the resource more strongly; Fig XXA). 
We find that consumer targeting strategies for different resource groups orient as spines emerging from a central point given by tau=0, i.e. a consumer that always targets the closest resource (black point).
The spines that branch out from the central cluster are those strategies targeting different resource groups.
Specialists on these resources are those farthest from the center, thus representing the most divergent strategies in this strategy-niche space.

How do different targeting strategies relate to estimates of consumer fitness, and is this predictive of consumer fitness in natural systems?
We next demonstrate how a simple measure of fitness for simulated consumers can - in principle - be used as an expectation for those observed at the Sevilleta.
Here we use the coefficient of variation (CV) of nitrogenous returns as a measure of fitness, such that lower values reflect smaller fluctuations relative to the mean (higher fitness), and higher values denote larger fluctuations relative to the mean (lower fitness).
Along the spines radiating from the central cluster, fitness values are generally consistent.
Across spines, this fitness landscape is roughly partitioned into those strategies resulting in higher fitness (lower CV; r) 
While the central spine (I) results in the lowest fitness (higher CV), those towards the central cluster (II) and along axes furthest from the central spine (III) demonstrate higher fitness (lower CV).

The diet of a hypothetical empirically-measured consumer is shown near the topmost spine (denoted by the `+').
Our establishment of the strategy-niche manifold in this context would thus offer an interpretation of this hypothetical consumer.
We would first note that the consumer is targeting C4 annual forbs, but that its strategy is one where roughly 50% of its foraging effort is oriented towards opportunistic resources.
Moreover, our foraging model predicts a fitness gradient with respect to this resource group, where increased specialization is expected to increase fitness, whereas increased reliance on opportunistic foraging decreases fitness (arrows).
Integrating individual-level ontogenetic information (Section xx) will allow us to directly assess the accuracy of this prediction.


## Confronting more realistic models with data (things to do)

While the simple foraging model that we describe here provides a useful heuristic for describing our approach, we additionally intend to incorporate modelling frameworks that include more realistic physical and biological constraints, which can be either combined to form a larger strategy-niche manifold or assessed independently.
For example, we intend to adapt a consumer foraging model described in Yeakel et al. (2020), where we use fitness maximization principles and stochastic dynamic programming to establish consumer strategies in environments with seasonal uncertainty.
In this case, our framework directly incorporates consumer energetic constraints, caching behaviors, and state-dependent foraging strategies, from which fitness is directly estimated.
Using more complex models such as this, in tandem with the simpler mechanics of the model described here, will not only enlarge our perspective of the universe of potential foraging strategies, but also enable us to pinpoint which biological and/or physical constraints play larger or smaller roles in driving consumer behaviors.

