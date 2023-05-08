Download Link: https://assignmentchef.com/product/solved-scicomp-project3
<br>
<h2>1.1         Importance of Rapid Convergence in Root Finding and Optimization</h2>

Finding a minimum energy configuration for a system of mutiple atoms is a widely applicaple for designing molecular systems in solid state physics, chemistry, medicine and other fields. Many applications will require tens to hundreds of atoms in the model. Calculating the energy involves solving the many-body Schrødinger equation, in which all the electrons interact with each other.

With hundreds of atoms and thousands of electrons, solving the many-body Schr¨odinger equation (even approximately) can take weeks of computetime to calculate a single energy point, whereby a full optimization can easily take many years of CPU-time to run. Even on a parallel computer with, say, 1000 CPU cores, optimizing the geometry of a single molecule can easily take months in “human time”.

When every evaluation of the energy function takes, for example, a CPU-week to compute (burning large amounts of CO<sub>2</sub>), the number of function evaluations matter: a problem that may be extra bad in quantum chemistry, but is common to many fields.

In this assignment, you will try your hand at optimizing a simple system of interacting atoms. However, we do not have months to sit and wait for computations to complete, so we will use an extremely simple approximation to the potential energy, the <em>Lennard-Jones</em>-potential. This approximation is appropriate for noble gas atoms, yielding a crude picture of the formation; yet we will pretend that it is a full <em>ab initio </em>quantum chemical energy calculation, and assume that it takes a week to complete. Thus, you will be asked to <em>count energy function evaluations </em>and report the time spent.

<h2>1.2         Minimizing the potential in a cloud of Argon atoms using the LennardJones potential</h2>

The Lennard-Jones potential works for systems of neutral atoms or molecules. It assumes that no new electronic bonds nor new molecules are formed, and that the potential can be modeled only with a short-range repulsion force (arising from the Pauli exclusion principle) and a long-range van der Waals attraction. This yields a classical description of the system, in which the neutral atoms or molecules move as classical particles, interacting only through this simple potential.

The Lennard-Jones potential can be written in several ways, but the most common is:

!

(1)

<em>v<sub>LJ</sub></em>(<em>r<sub>ij</sub></em>) is the potential at distance <em>r<sub>ij </sub></em>= k<strong>x</strong><em><sub>i </sub></em>− <strong>x</strong><em><sub>j</sub></em>k<sub>2 </sub>between two atoms <em>i </em>and <em>j</em>. The repulsive force is stronger than the van der Waals-force, but decreases more rapidly with distance. Thus, for a pair of atoms, there is a distance where the total potential is minimal, the potential well. <em> </em>is this minimal potential between two atoms, and <em>σ </em>is the inter-particle distance where the potential is zero: <em>v<sub>LJ</sub></em>(<em>σ</em>) = 0. The constants <em> </em>and <em>σ </em>are generally found by fitting to proper ab initio calculations, or using experimental data.

The full potential LJ-energy of <em>N </em>neutral particles is the sum of pair-potentials:

<em>N           N</em>

<table width="564">

 <tbody>

  <tr>

   <td width="547"><em>V</em><em>LJ</em>(<strong>X</strong>) = X X <em>v</em><em>LJ</em>(<em>r</em><em>ij</em>)<em>i</em>=1 <em>j</em>=<em>i</em>+1where <strong>X </strong>∈ R<em><sup>N</sup></em><sup>×3 </sup>is the matrix of coordinates for the <em>N </em>particles.In this assignment, we will be modeling argon, for which</td>

   <td width="17">(2)</td>

  </tr>

 </tbody>

</table>

<em>σ </em>= 3<em>.</em>401˚A and  997kJ/mol                                                         (3)

<h2>1.3         Array Programming</h2>

Try to use array programming when programming with Python/Numpy or Matlab, as the speed difference can be up to a factor of hundreds. As an example of how to use array programming, here is how a distance matrix (to be used in the LJ-potential) can be calculated with array operations using NA = newaxis (as I explain further <a href="http://www.nbi.dk/~avery/teaching/scicomp/numpy-programming.pdf">here)</a>:

<em>D</em><em>ij </em>= k<strong>x</strong><em>i </em>− <strong>x</strong><em>j</em>k2 = k<strong>x</strong><em>i</em><strong>j </strong>− <strong>x</strong><strong>i</strong><em>j</em>k2                                                                                                 (4)

In Eq. (4), <strong>x</strong><em><sub>i </sub></em>is extended with a newaxis <strong>j </strong>and <strong>x</strong><em><sub>j </sub></em>with a newaxis <strong>i </strong>to put them on the same footing, so that the pair differences can be calculated as elementwise subtraction <strong>x</strong><em><sub>i</sub></em><strong><sub>j </sub></strong>− <strong>x</strong><strong><sub>i</sub></strong><em><sub>j</sub></em>. In Numpy code, it looks as follows:

<em># points:                                   (N,3)</em>−<em>array of (x,y,z) coordinates for N points</em>

<em># distance(points): returns (N,N)</em>−<em>array of inter</em>−<em>point distances:</em>

def distance( points ):

displacement = points[:,NA] − points[NA,:] return sqrt( sum(displacement∗displacement, axis=−1) )

<h1>2           Questions for Week 4: Solving Nonlinear Equations</h1>

<h2>2.1         Solving Nonlinear Equations in 1D</h2>

<ol>

 <li>Write a function that computes the Lennard-Jones energy for a collection of particles described by an (<em>N,</em>3)-array points. The resulting function should not take more arguments than this, as it will be called with by your root-finding functions below. In Python, you can do this using a programming construct called <em>closures </em>as follows:</li>

</ol>

def LJ(sigma,epsilon): def V( points ):

…implementation that depends on sigma and epsilon… return V

Then) will return a function specialized to <em>σ </em>and <em> </em>such that <em>V </em>(points) returns the total Lennard-Jones energy for the system.

Demonstrate your solution by making 1) a plot of the potential between two Ar atoms, one placed at <strong>x</strong><sub>0 </sub>= (<em>x,</em>0<em>,</em>0) the other at <strong>x</strong><sub>1 </sub>= (0<em>,</em>0<em>,</em>0) with <em>x </em>ranging from 3 to 11; and 2) a plot of the LJ-potential in the same range with two extra points added: <strong>x</strong><sub>2 </sub>= (14<em>,</em>0<em>,</em>0) and <strong>x</strong><sub>3 </sub>= (7<em>,</em>3<em>.</em>2<em>,</em>0).

<ol>

 <li>Write a bisection root finding function x, ncalls = bisectionroot(f,a,b,tolerance=1e−13) that finds x such that <em>f</em>(x) = 0 given a bracket <em>x </em>∈ [<em>a</em>;<em>b</em>], and counts the number of calls to the function <em>f</em>. (A reasonable length is 8-12 lines of code).</li>

</ol>

In this assignment, let the convergence test be on how close we get <em>f</em>(<em>x</em>) to zero. Test it to find the zero of the the LJ-potential between two argon atoms as a function of interatomic distance, and verify that you get <em>x </em>= <em>σ</em>. How many calls to the energy function were needed to get from the start bracket [<em>a,b</em>] = [2<em>,</em>6] to |<em>f</em>(<em>x</em>)| <em>&lt; </em>10<sup>−13</sup>?

<ol>

 <li>The derivative of the pair-potential) is</li>

</ol>

Write a Newton-Rhapson solver x, ncalls = newtonroot(f,df,x0,tolerance,maxiterations), (a reasonable length is 4-8 lines of code), and test it in the same way as above. For simplicity, assume that a call to the derivative df has the same cost as a call to f. How many calls were needed to get from <em>x</em><sub>0 </sub>= 2 to |<em>f</em>(<em>x</em><sup>∗</sup>)| <em>&lt; </em>10<sup>−12</sup>, i.e., 12 decimal digits for <em>x</em><sup>∗ </sup>after the comma?

<ol>

 <li>Make a combination of Newton-Rhapson and bisection that is <em>guaranteed to converge</em>, but takes advantage of the quadratic convergence of Newton-Rhapson iteration, and test it on the same example. How many calls to the LJ-energy function was needed to get from <em>x</em><sub>0 </sub>= 2, [<em>a,b</em>] = [2<em>,</em>6] to obtain |<em>f</em>(<em>x</em><sup>∗</sup>)| <em>&lt; </em>10<sup>−13</sup>?</li>

</ol>

<strong>Note: </strong><em>If you have trouble completing this step, simply skip it and move on to the remaining questions, which only requires your bisection root solver to work. You can always return to solve it once you have completed tasks (e) and (f).</em>

<h2>2.2         Solving <em>N</em>-dimensional Nonlinear Equations</h2>

Using the chain rule for derivatives <em>f</em>(<em>g</em>(<em>x</em>))<sup>0 </sup>= <em>f</em><sup>0</sup>(<em>g</em>(<em>x</em>))<em>g</em><sup>0</sup>(<em>x</em>), we can write down the derivative of a pair-potential with respect to the position of one of the particles:<sup>1</sup>

(5)

This we can then use to find an expression for the gradient of the total LJ-energy:

(6)

The position of all the particles is a point <strong>X </strong>∈ R<sup>3<em>N</em></sup>, which we organize as a (<em>N,</em>3)-array, so that <strong>x</strong><em><sub>i </sub></em>= (<em>x<sub>i</sub>,y<sub>i</sub>,z<sub>i</sub></em>). The total potential energy is a function <em>V<sub>LJ </sub></em>: R<sup>3<em>N </em></sup>→ R, and the gradient taken at any particular configuration of particle positions is hence a 3<em>N</em>-dimensional vector: <em>V<sub>LJ</sub></em>(<strong>X</strong>) ∈ R<sup>3<em>N</em></sup>. Thus, the gradient to <em>V<sub>LJ </sub></em>is a function ∇<em>V<sub>LJ </sub></em>: R<sup>3<em>N </em></sup>→ R<sup>3<em>N</em></sup>. The negative of the gradient is the <em>force </em>acting on the system, and its direction is that in which the potential decreases most rapidly. But notice that it is a 3<em>N</em>-dimensional direction: it acts on <em>all </em>the particles at once.

We can write it down in Python as <a href="http://www.nbi.dk/~avery/teaching/scicomp/2020/LJhelperfunctions.py">follows</a> (in Matlab, it would look very similar):

def LJgradient(sigma, epsilon): def gradV(X): d = X[:,NA] − X[NA,:] <em># (N,N,3) displacement vectors </em>r = sqrt( sum(d∗d,axis=−1) ) <em># (N,N) distances </em>fill diagonal(r,1) <em># Don’t divide by zero</em>

T = 6∗(sigma∗∗6) ∗ (r∗∗−7) − 12∗(sigma∗∗12) ∗ (r∗∗−13) <em># (N,N)</em>−<em>matrix of r</em>−<em>derivatives</em>

<em># Using the chain rule, we turn the (N,N)</em>−<em>matrix of r</em>−<em>derivatives into # the (N,3)</em>−<em>array of derivatives to Cartesian coordinate: the gradient.</em>

<em># (Automatically sets diagonal to (0,0,0) = X[i]</em>−<em>X[i]) </em>u = d/r[:,:,NA]         <em># u is (N,N,3)</em>−<em>array of unit vectors in direction of X[i] </em>− <em>X[j] </em>return 4∗epsilon∗sum(T[:,:,NA] ∗ u, axis=1)

return gradV

<sup>1</sup>NB: You don’t need to understand these derivations to solve the problem, as the actual code for the gradient is provided to you. The derivations are here to show you how to do it yourself in the future.

Finding the minima in R<sup>3<em>N </em></sup>of the potential involves finding points where its 3<em>N</em>-dimensional gradient is <strong>0</strong>. An important component needed in next week’s work is called <em>line search</em>, where we search for a point along a single direction <strong>d </strong>∈ R<sup>3<em>N </em></sup>at which the gradient <em>along this line </em>is zero. That is, we want to start in a point <strong>X</strong><sub>0 </sub>∈ R<sup>3<em>N</em></sup>, and then find <em>α </em>∈ [0;<em>b</em>] such that 0 = <strong>d </strong>· ∇<em>V<sub>LJ</sub></em>(<strong>X</strong><sub>0 </sub>+ <em>α</em><strong>d</strong>), i.e., the gradient has no component in the direction of <strong>d</strong>. This finds an optimum for the one-dimensional function <em>V<sub>LJ</sub></em>(<strong>X</strong><sub>0 </sub>+ <em>α</em><strong>d</strong>).

<ol>

 <li>Look at the gradient of the 2-particle system in question (a) with <strong>x</strong><sub>1 </sub>= (0<em>,</em>0<em>,</em>0) and <strong>x</strong><sub>0 </sub>= (<em>x,</em>0<em>,</em>0), <em>x </em>∈ [3;10]. Why are exactly two components nonzero? Why are they equal and opposite? Plot the nonzero component for the derivative of the <em>x</em>-coordinate of <strong>x</strong><sub>0 </sub>(0<em>,</em>0coordinate of gradient) together with the potential, and notice the relationship between the zero of the derivative and the minimum of the potential.</li>

</ol>

Next look at the gradient for the 4-particle system from (a) at one of the minima of your plot. Why is the gradient not zero?

<ol>

 <li>Write a function x0, ncalls = linesearch(F,X0, d, alphamax, tolerance, maxiterations) that takes a function <strong>F </strong>: R<em><sup>N</sup></em><sup>×3 </sup>→ R<em><sup>N</sup></em><sup>×3</sup>, a start position <strong>X</strong><sub>0 </sub>∈ R<em><sup>N</sup></em><sup>×3 </sup>and finds the zero along the line-segment <strong>X</strong><sub>0 </sub>+ <em>α</em><strong>d </strong>of <strong>d </strong> <strong>F</strong>(<strong>X</strong><sub>0 </sub>+ <em>α</em><strong>d</strong>).<a href="#_ftn1" name="_ftnref1"><sup>[1]</sup></a> Use your bisection solver at this point, as using Newton-Rhapson requires second derivatives when <strong>F </strong>is the gradient.</li>

</ol>

Test your function by finding the mimimum along <strong>X</strong><sub>0 </sub>+ <em>α</em><strong>d </strong>with

<strong>X</strong>

and <strong>d </strong>= −∇<em>V<sub>LJ</sub></em>(<strong>X</strong><sub>0</sub>), and with <em>α </em>∈ [0;1].

<a href="#_ftnref1" name="_ftn1">[1]</a> You can use either a closure or a lambda expression to define the one-dimensional restriction of <strong>F </strong>to the line.