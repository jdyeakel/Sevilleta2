# Diffusion Map of the Dietary Niche

## From N-D Hamiltonian Space to Niche Manifold

We will now address how we will we compare the foraging behavior of organisms operating within complex seasonal environments, and relate these behaviors to measures of fitness in a framework that is investigative and potentially predictive.
Hutchinson famously proposed an N-dimensional hypervolume to describe an organism's niche, however we must suppose that most of this potential volume cannot be realized due to physical, behavioral, or environmental constraints. 
We can therefore imagine a process by which the potential space can be observed apart from that which is nonphysical or unrealizable for a particular suite of organisms to gain a clearer understanding of species' roles within a larger community.
In the context of the dietary niche, mechanistic foraging models provide a means by which to quantify the potential niche space given a range of constraints, as the consequences of a large range of foraging strategies can be systematically assessed -- in the albeit simplified manner of mathematical models.
If the constraints captured by the foraging model play a central role in shaping the foraging behaviors of consumers in natural systems, the niche-space enumerated by the model can be used as a lens by which the larger community can be understood.
However, even a simple foraging model can result in an enormous complexity of foraging strategies, and a central challenge in understanding both model and natural ecological systems resides in excising low-dimensional relationships from high-dimensional data, particularly when such data are generated from interactions with nonlinear dependencies.

Here we propose to devise a series of foraging models that range from simple to complex, to explore and enumerate the dietary consequences of a large range of foraging strategies.
We will then use this range of strategies to construct a low-dimensional realization of a high-dimensional niche space using diffusion mapping techniques.
Diffusion mapping techniques constitute a class of data analytics capable of reconstructing nonlinear relationships governing high-dimensional datasets when the generative processes underlying the data is unknown.
We will use this perspective to reconstruct a low-dimensional embedding of consumer strategies, or *niche manifold*, capturing the structural similarities of a range of foraging strategies.
Understanding how different strategies relate to one another when the mechanics are simple and well-defined provides a null expectation by which to classify and evaluate the more complex strategies observed among consumers in natural systems.
Our procedure thus consists of three parts: to 1) build a class of foraging models to systematically capture a range of foraging strategies given known physical and biological constraints; 2) reconstruct a low-dimensional embedding of the range of modeled foraging behaviors to capture the associated *niche manifold*; 3) measure the diets of local consumers and assess their ecological roles by their proximity to modeled strategies on the niche manifold; 4) explore predictive...

Because both the underlying generative models and resulting foraging strategies of model systems within the embedding are known, in addition to assessment, we can use the constructed niche manifold in a predictive capacity.
For instance, the foraging constraints that give rise to different strategies also allow us to assess the fitness consequences of those strategies.
By relating the fitness consequences of simulated strategies occupying different locations along the niche manifold, we can associate empirically measured strategies with modeled values nearby.
The resultant fitness landscape assessed along the niche manifold can thus be used to evaluate and ultimately predict life history characteristics of consumers foraging in natural systems.
Thus, by incorporating foraging data and estimates of fitness from Sevilleta consumers alongside those simulated within the context of the niche manifold, we will directly assess the predictive value of the manifold niche concept in a natural community.

## A simple foraging null model

We present a minimal consumer foraging model to illustrate how our framework can be used to uncover a simulated consumer's niche manifold -- which we will identify using diffusion mapping techniques -- and treat this as a null expectation by which we will subsequently evaluate empirical consumer strategies.
We note that mechanistic foraging models of arbitrary complexity could be used to build this null expectation, and while we will identify in what direction we aim to build additional complexity, a simpler dynamic best suits an illustration of the core approach.

We simulate a consumer foraging in a seasonal environment given a set of resource functional groups (see xx) with differing spatial distributions during summer and winter seasons, such that mu_si is the mean encounter rate where s denotes season and i=1...N functional groups.
A consumer of mass M forages within this landscape, targeting a particular functional group with weight tau.
The targeting weight tau means that, for each consumer-resource interaction, the consumer will find and acquire its targeted resource group with probability tau, and target the *closest* resource group with probably 1-tau, regardless of preference.
Once a consumer-resource interaction is drawn, the consumer travels the distance to the resource with velocity v(M), and assimilates both bulk energy (kJ) and nitrogen content based on the resource's energy density and nitrogen concentration, respectively.
Consumer-resource interactions continue until a predetermined time threshold for the foraging bout, tmax, is reached, whereupon 

## The niche manifold of the null model

## Confronting the model with data (things to do)