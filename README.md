Download Link: https://assignmentchef.com/product/solved-fys4150-project-2
<br>



<h1>Eigenvalue problems, from the equations of a buckling beam to Schroedinger’s equation for two electrons in a threedimensional harmonic oscillator well</h1>

<strong>Introduction. </strong>The aim of this project is to develop your own code for solving eigenvalue problems. The matrix to diagonalize is the same as the one we encountered in project one, a so-called tridiagonal Toeplitz matrix. This matrix has analytical eigenpairs (eigenvalues and eigenvectors) and gives us an excellent testing ground for our algorithms. In this project we will develop an eigenvalue solver based on Jacobi’s method. The project will also introduce you to units tests and we will compare our results against other eigenvalue solvers (from LAPACK and/or numpy).

This project aims also at introducing to you the concept of <a href="https://www.springer.com/us/book/9783319327259">scaling of equations</a><a href="https://www.springer.com/us/book/9783319327259">. </a>This means often either making various variables dimensionless or introducing units which are more convenient.

We will start with the two-point boundary value problem of a buckling beam or a spring fastened at both ends. This is one of the problems which has analytical solutions. Thereafter, by simply adding a new variable along the diagonal elements, we can study quantum mechanical problems. In particular, we will study the harmonic oscillator problem in three dimensions, with one or two electrons. For the latter case we can study the role of the repulsive Coulomb interaction and extract interesting physics results. For selected frequencies, this interacting two-electron problem exhibits analytical solutions, one of the few cases of an interacting system where wecan find analytical solutions. See <a href="https://prola.aps.org/abstract/PRA/v48/i5/p3561_1">M. </a><a href="https://prola.aps.org/abstract/PRA/v48/i5/p3561_1">Taut, Phys. Rev. A 48, 3561 (1993)</a> for the derivation of analytical expressions for the eigenpairs..

Electrons confined in small areas in semiconductors, so-called quantum dots, form a hot research area in modern solid-state physics, with applications spanning from such diverse fields as quantum nano-medicine to the contruction of quantum gates. The article on quantum computing with quantum dots by <a href="https://journals.aps.org/pra/abstract/10.1103/PhysRevA.57.120">Loss and DiVincenzo</a> is an excellent read for those interested in this exciting topic.

<em> </em>c 1999-2019, “Computational Physics I

FYS3150/FYS4150″:”http://www.uio.no/studier/emner/matnat/fys/FYS3150/indexeng.html”. Released under CC Attribution-NonCommercial 4.0

license

<strong>The buckling beam problem, a classical wave function problem in one dimension. </strong>We start with the following differential equation, namely

<em>d</em><sup>2</sup><em>u</em>(<em>x</em>)

<em>γ    </em> = <em>−Fu</em>(<em>x</em>)<em>, dx</em><sup>2</sup>

where <em>u</em>(<em>x</em>) is the vertical displacement of the beam in the <em>y </em>direction. The beam has length <em>L</em>, <em>x ∈ </em>[0<em>,L</em>] and <em>F </em>is a force applied at (<em>L,</em>0) in the direction towards the origin. The parameter <em>γ </em>is a constant defined by properties like the rigidity of the beam. We apply again so-called Dirichlet boundary conditions and set <em>u</em>(0) = <em>u</em>(<em>L</em>) = 0.

In this specific case two of the parameters <em>γ</em>, <em>F </em>and <em>L </em>are known. As an example, assume we know <em>F </em>and <em>L</em>. Then the eigenvalue problem we set up below will allow us to find <em>γ</em>.

We define a dimensional variable

<em>x ρ </em>=        <em>,</em>

<em>L</em>

meaning that we have <em>ρ ∈ </em>[0<em>,</em>1]. By reordering the equation as

<em>d</em><sup>2</sup><em>u</em>(<em>ρ</em>)            <em>FL</em><sup>2</sup>

= <em>− u</em>(<em>ρ</em>) = <em>−λu</em>(<em>ρ</em>)<em>, dρ</em><sup>2                    </sup><em>γ</em>

with <em>λ </em>= <em>FL</em><sup>2</sup><em>/γ </em>we have an equation that when discretized, becomes an eigenvalue problem. We use the same expression for the second derivative of a function <em>u </em>as we did in project 1, namely

<em><sub>00                </sub>u</em>(<em>ρ </em>+ <em>h</em>) <em>− </em>2<em>u</em>(<em>ρ</em>) + <em>u</em>(<em>ρ − h</em>)                  <sub>2</sub>

<em>u </em>=                               + <em>O</em>(<em>h </em>)<em>,     </em>(1) <em>h</em>2

where <em>h </em>is our step. Next we define minimum and maximum values for the variable <em>ρ</em>, <em>ρ</em><sub>min </sub>= 0 and <em>ρ</em><sub>max </sub>= 1, respectively. With a given number of mesh points, <em>N</em>, we define the step length <em>h </em>as, with <em>ρ</em><sub>min </sub>= <em>ρ</em><sub>0 </sub>and <em>ρ</em><sub>max </sub>= <em>ρ<sub>N</sub></em>,

<em>ρ</em><em>N − ρ</em>0 <em>h </em>=     <em>.</em>

<em>N</em>

The value of <em>ρ </em>at a point <em>i </em>is then

<em>ρ<sub>i </sub></em>= <em>ρ</em><sub>0 </sub>+ <em>ih                    i </em>= 1<em>,</em>2<em>,…,N.</em>

We can rewrite the differential equation for a value <em>ρ<sub>i </sub></em>as

<em>−u</em>(<em>ρ</em><em>i </em>+ <em>h</em>) <em>− </em>2<em>u</em>(2<em>ρ</em><em>i</em>) + <em>u</em>(<em>ρ</em><em>i − h</em>) = <em>λu</em>(<em>ρ<sub>i</sub></em>)<em>,</em>

<em>h</em>

or in a more compact way as

<em>−u</em><em>i</em>+1 <em>− </em>2<em>u</em>2<em>i </em>+ <em>u</em><em>i−</em>1 <sub>= <em>λu</em><em>i</em></sub><em><sub>.</sub></em>

<em>h</em>

Following our approach from project 1, we can rewrite this equation in a more a general form, but now as an eigenvalue problem, that is

<sup> </sup><em>d        a       </em>0       0       <em>…       </em>0            0 <sup> </sup><em>u</em><sub>1 </sub><sup>                </sup> <em>u</em><sub>1 </sub>

<sub> </sub><em>a        d       a       </em>0       <em>…       </em>0            0 <sub> </sub><em>u</em><sub>2 </sub><sub>                </sub> <em>u</em><sub>2 </sub>

 0        <em>a       d       a       </em>0      <em>…      </em>0  <em>u</em>3  <em>λ</em> <em>…u</em><sup>3 </sup><sup></sup><em>.              </em>(2)



<sup></sup><em>…     …      …      …      …      …        …</em><sup></sup><sub></sub><sup></sup><sub> </sub><em>… </em><sup></sup><sub> </sub><sup>=</sup>

 0       <em>…      …      …       a       d       a </em><em>u</em><em>N−</em>2         <em>u</em><em>N−</em>2



0      <em>…      …      …      …       a       d          u<sub>N−</sub></em><sub>1                            </sub><em>u</em><em>N−</em>1

As in project 1, we have not included the endpoints <em>u</em><sub>0 </sub>and <em>u<sub>N</sub></em>. We have defined <em>d </em>= 2<em>/h</em><sup>2 </sup>and the non-diagonal ones as <em>a </em>= <em>−</em>1<em>/h</em><sup>2</sup>. This eigenvalue problem has analytical eigenpairs, with eigenvalues given as

<em>jπ λ<sub>j </sub></em>= <em>d </em>+ 2<em>a</em>cos( ) <em>j </em>= 1<em>,</em>2<em>,…N. N </em>+ 1

<strong>Project 2 a): Mathematical intermezzo. </strong>A unitary transformation preserves the orthogonality of the obtained eigenvectors. To see this consider first a basis of vectors <strong>v</strong><em><sub>i</sub></em>,

<em>v</em><em>i</em>1

<strong>v</strong><em>i </em>= <sub></sub><em>……</em><sub></sub>

<em>v</em><em>in</em>

We assume that the basis is orthogonal, that is

<strong>v</strong><em>jT</em><strong>v</strong><em>i </em>= <em>δ</em><em>ij.</em>

Show that an orthogonal or unitary transformation

<strong>w</strong><em><sub>i </sub></em>= <strong>Uv</strong><em><sub>i</sub>,</em>

preserves the dot product and orthogonality.

<strong>Project 2 b): Setting up a code for tridiagonal Toeplitz matrix. </strong>Your task now is to write a function which implements Jacobi’s rotation algorithm (see Lecture notes chapter 7) in order to solve Eq. (2). However, the first simple check is to set up the matrix to diagonalize for a given <em>N </em>and use either armadillo’s or numpy’s functions for diagonalizing a matrix. You can then check that you obtain the analytical eigenvalues.

For Jacobi’s method, we define the quantities tan<em>θ </em>= <em>t </em>= <em>s/c</em>, with <em>s </em>= sin<em>θ </em>and <em>c </em>= cos<em>θ </em>and

<em>a</em><em>ll − a</em><em>kk </em>cot2<em>θ </em>= <em>τ </em>= <em>.</em>

2<em>a<sub>kl</sub></em>

We can then define the angle <em>θ </em>so that the non-diagonal matrix elements of the transformed matrix <em>a<sub>kl </sub></em>become non-zero and we obtain the quadratic equation

(using cot2<em>θ </em>= 1<em>/</em>2(cot<em>θ − </em>tan<em>θ</em>)

<em>t</em><sup>2 </sup>+ 2<em>τt − </em>1 = 0<em>,</em>

resulting in

<em>t </em>= <em>−τ ± </em><sup>p</sup>1 + <em>τ</em><sup>2</sup><em>,</em>

and <em>c </em>and <em>s </em>are easily obtained via

<em>c </em><em>,</em>

and <em>s </em>= <em>tc</em>.

How many similarity transformations are needed before you reach a result where all non-diagonal matrix elements are essentially zero? Try to estimate the number of transformations and extract a behavior as function of the dimensionality of the matrix. Compare your results with the analytical ones.

You can check your results against the Armadillo function for solving eigenvalue problems. The armadillo function <em>eigsys </em>can be used to find eigenvalues and eigenvectors. A Python program which solves this part of the project is available under the <a href="https://compphysics.github.io/ComputationalPhysics/doc/pub/projectwriting/html/projectwriting.html">project writing slides</a><a href="https://compphysics.github.io/ComputationalPhysics/doc/pub/projectwriting/html/projectwriting.html">.</a>

Comment your results (here you could for example compute the time needed for both algorithms for a given dimensionality of the matrix).

<strong>Project 2 c): Implementing tests in your code. </strong>In this project (and later ones as well) we will implement so-called <strong>unit </strong>tests. Our unit tests are mainly meant to test mathematical properties of our algorithm. During the development phase of a program it is useful to devise tests that your program should pass. One of these is to make sure that for a given simple test matrix (say a 5 <em>× </em>5 matrix) our algorithm for searching for the largest non-diagonal element always returns the correct answer. Furthermore, for a known simple matrix, irrespective of changes made, we should always get the same and correct eigenvalues. Another test is to make sure that the orthogonality shown in exercise (a) is preserved. Try to figure out other tests your code should pass, based either on the mathematical properties of the algorithms or more program specific tests. Implement at least two unit tests in this project.

<strong>Extending our machinery to quantum mechanics.      </strong>Here we will assume that these electrons move in a three-dimensional harmonic oscillator potential (they are confined by for example quadrupole fields) and repel each other via the static Coulomb interaction. We assume spherical symmetry. You don’t need to think of quantum mechanics when solving this project, look at it as another diagonalization problem.

We are first interested in the solution of the radial part of Schroedinger’s equation for one electron. This equation reads

<em>−</em><em><sup>d </sup>r</em><sup>2 </sup><em><sup>d </sup>− </em><em><sup>l</sup></em><sup>(<em>l </em>+ 1)</sup><em>R</em>(<em>r</em>) + <em>V </em>(<em>r</em>)<em>R</em>(<em>r</em>) = <em>ER</em>(<em>r</em>)<em>.</em>

2<em>m      r</em><sup>2 </sup><em>dr       dr            r</em><sup>2</sup>

In our case <em>V </em>(<em>r</em>) is the harmonic oscillator potential (1<em>/</em>2)<em>kr</em><sup>2 </sup>with <em>k </em>= <em>mω</em><sup>2 </sup>and <em>E </em>is the energy of the harmonic oscillator in three dimensions. The oscillator frequency is <em>ω </em>and the energies are

3

<em>E<sub>nl </sub></em><em>     ,</em>

with <em>n </em>= 0<em>,</em>1<em>,</em>2<em>,… </em>and <em>l </em>= 0<em>,</em>1<em>,</em>2<em>,…</em>.

Since we have made a transformation to spherical coordinates it means that <em>r ∈ </em>[0<em>,∞</em>). The quantum number <em>l </em>is the orbital momentum of the electron. Then we substitute <em>R</em>(<em>r</em>) = (1<em>/r</em>)<em>u</em>(<em>r</em>) and obtain

~2 <em>d</em>2                                                 <em>l</em>(<em>l </em>+ 1) ~2

<em>−  </em><em>u</em>(<em>r</em>) +       <em>V </em>(<em>r</em>) +    <em>      u</em>(<em>r</em>) = <em>Eu</em>(<em>r</em>)<em>.</em>

2<em>m dr</em><sup>2                                                          </sup><em>r</em><sup>2           </sup>2<em>m</em>

The boundary conditions are <em>u</em>(0) = 0 and <em>u</em>(<em>∞</em>) = 0.

We introduce a dimensionless variable <em>ρ </em>= (1<em>/α</em>)<em>r </em>where <em>α </em>is a constant with dimension length and get

~2         <em>d</em>2                                                  <em>l</em>(<em>l </em>+ 1) ~2

<em>−  </em><em>u</em>(<em>ρ</em>) +      <em>V </em>(<em>ρ</em>) +   <em>      u</em>(<em>ρ</em>) = <em>Eu</em>(<em>ρ</em>)<em>.</em>

2<em>mα</em><sup>2 </sup><em>dρ</em><sup>2                                                           </sup><em>ρ</em><sup>2          </sup>2<em>mα</em><sup>2</sup>

We will set in this project <em>l </em>= 0. Inserting <em>V </em>(<em>ρ</em>) = (1<em>/</em>2)<em>kα</em><sup>2</sup><em>ρ</em><sup>2 </sup>we end up with

<em>− </em>~2 <em>d</em>2 <em>u</em>(<em>ρ</em>) + <em>kα</em><sup>2</sup><em>ρ</em><sup>2</sup><em>u</em>(<em>ρ</em>) = <em>Eu</em>(<em>ρ</em>)<em>.</em>

2<em>mα</em><sup>2 </sup><em>dρ</em><sup>2                        </sup>2

We multiply thereafter with 2<em>mα</em><sup>2</sup><em>/</em>~<sup>2 </sup>on both sides and obtain

<em>− </em><em>d</em>2 <em>u</em>(<em>ρ</em>) + <em>mk</em><em>α</em>4<em>ρ</em>2<em>u</em>(<em>ρ</em>) = 2<em>mα</em>2 <em>Eu</em>(<em>ρ</em>)<em>. dρ</em>2 ~2 ~2

The constant <em>α </em>can now be fixed so that

<em>mk </em><sub>4 </sub><em>α </em>= 1<em>,</em>

~2

or

<em>α </em>=    <em>         . mk</em>

Defining

<em>E,</em>

we can rewrite Schroedinger’s equation as

<em>− </em><em>d</em>2 <em>u</em>(<em>ρ</em>) + <em>ρ</em><sub>2</sub><em>u</em>(<em>ρ</em>) = <em>λu</em>(<em>ρ</em>)<em>. dρ</em><sup>2</sup>

This is the first equation to solve numerically. In three dimensions the eigenvalues for <em>l </em>= 0 are <em>λ</em><sub>0 </sub>= 3<em>,λ</em><sub>1 </sub>= 7<em>,λ</em><sub>2 </sub>= 11<em>,….</em>

We define minimum and maximum values for the variable <em>ρ</em>, <em>ρ</em><sub>min </sub>= 0 and <em>ρ</em><sub>max</sub>, respectively. You need to check your results for the energies against different values <em>ρ</em><sub>max</sub>, since we cannot set <em>ρ</em><sub>max </sub>= <em>∞</em>.

With a given number of mesh points, <em>N</em>, we define the step length <em>h </em>as, with <em>ρ</em>min = <em>ρ</em>0 and <em>ρ</em>max = <em>ρ</em><em>N</em>,

<em>ρ</em><em>N − ρ</em>0 <em>h </em>=     <em>.</em>

<em>N</em>

The value of <em>ρ </em>at a point <em>i </em>is then

<em>ρ<sub>i </sub></em>= <em>ρ</em><sub>0 </sub>+ <em>ih                    i </em>= 1<em>,</em>2<em>,…,N.</em>

We can rewrite the Schroedinger equation for a value <em>ρ<sub>i </sub></em>as

<em>−</em><em>u</em>(<em>ρ</em><em>i </em>+ <em>h</em>) <em>− </em>2<em>u</em>(2<em>ρ</em><em>i</em>) + <em>u</em>(<em>ρ</em><em>i − h</em>) + <em>ρ</em><sub>2<em>i</em></sub><em>u</em>(<em>ρ<sub>i</sub></em>) = <em>λu</em>(<em>ρ<sub>i</sub></em>)<em>, h</em>

or in a more compact way

<em>−u</em><em>i</em>+1 <em>− </em>2<em>u</em>2<em>i </em>+ <em>u</em><em>i−</em>1 + <em>ρ</em>2<em>i</em><em>u</em><em>i </em>= <em>−</em><em>u</em><em>i</em>+1 <em>− </em>2<em>u</em>2<em>i </em>+ <em>u</em><em>i−</em>1 + <em>V</em><em>iu</em><em>i </em>= <em>λu</em><em>i, h          h</em>

where <em>V<sub>i </sub></em>= <em>ρ</em><sup>2</sup><em><sub>i </sub></em>is the harmonic oscillator potential. We define first the diagonal matrix element

2

<em>d</em><em>i </em>= <em>h</em>2 + <em>V</em><em>i,</em>

and the non-diagonal matrix element

1

<em>e</em><em>i </em>= <em><sup>− </sup></em>2 <em>. </em><em><sub>h</sub></em>

In this case the non-diagonal matrix elements are given by a mere constant. <em>All non-diagonal matrix elements are equal</em>. With these definitions the Schroedinger

equation takes the following form

<em>d</em><em>iu</em><em>i </em>+ <em>e</em><em>iu</em><em>i−</em>1 + <em>e</em><em>iu</em><em>i</em>+1 = <em>λu</em><em>i,</em>

<table width="459">

 <tbody>

  <tr>

   <td colspan="8" width="459">where <em>u<sub>i </sub></em>is unknown. We can write the latter equation as a matrix eigenvalue</td>

  </tr>

  <tr>

   <td width="59">problem</td>

   <td width="29"> </td>

   <td width="29"> </td>

   <td width="29"> </td>

   <td width="61"> </td>

   <td width="44"> </td>

   <td width="192"> </td>

   <td width="17"> </td>

  </tr>

  <tr>

   <td width="59"><em>d</em><sub>1</sub><em>e</em><sub>1</sub> 0<sup></sup><em>…</em>  00</td>

   <td width="29"><em>e</em><sub>1</sub><em>d</em><sub>2 </sub><em>e</em><sub>2 </sub><em>… …</em><em>…</em></td>

   <td width="29">0 <em>e</em><sub>2</sub><em>d</em><sub>3</sub><em>… …</em><em>…</em></td>

   <td width="29">00 <em>e</em><sub>3</sub><em>… …</em><em>…</em></td>

   <td width="61"><em>…</em><em>… </em>0<em>…</em><em>…e</em><em>N−</em>3<em>…</em></td>

   <td width="44">0 0<em>…</em><em>…</em><em>d</em><em>N−</em>2 <em>e</em><em>N−</em>2</td>

   <td width="192">0  <em>u</em><sub>1 </sub>                            <em>u</em><sub>1 </sub>0  <em>u</em>2                             <em>u</em>2 0  <em>… </em> <em>λ</em> <em>…… </em><em>.</em><em>… </em> <em>… </em> = <em>e</em><em>N−</em>2 <em>… </em>                 <em>… </em><em>d</em><em>N−</em>1                 <em>u</em><em>N−</em>1                           <em>u</em><em>N−</em>1</td>

   <td width="17">(3)</td>

  </tr>

 </tbody>

</table>

Note here that we have not included the endpoints since we assume that they are knowm (as we did in project 1). The matrix is also symmetric in this case.

<strong>Project 2 d): Quantum dots in three dimensions, one electron. </strong>Add the harmonic oscillator potential to your tridiagonal Toeplitz matrix from 2a-2c and diagonalize the matrix. Study the results as functions of the number of integration points <em>N </em>and your approximation to <em>ρ</em><sub>max</sub>. The analytical results with our scaling for the one-electron energies are <em>λ </em>= 3<em>,</em>7<em>,</em>11<em>,</em>15<em>,…</em>. How many integration points do you need in order to reproduce the analytical results with say four leading digits after the decimal point?

You can reuse the code you wrote for part (b), but you need to add the potential <em>ρ</em><sup>2 </sup>to the diagonal elements.

<strong>Project 2 e): Quantum dots in three dimensions, two electrons. </strong>Note: if you are not interested in the explicit quantum physics case here, you can replace this part with either 2g or 2h below. The latter ones are optional exercises if you do this.

We will now study two electrons in a harmonic oscillator well which also interact via a repulsive Coulomb interaction. Let us start with the single-electron equation written as

~<sup>2 </sup><em>d</em><sup>2</sup>

<em>− </em><em>u kr</em><sup>2</sup><em>u</em>(<em>r</em>) = <em>E</em><sup>(1)</sup><em>u</em>(<em>r</em>)<em>, </em>2<em>m dr</em><sup>2 </sup>2

where <em>E</em><sup>(1) </sup>stands for the energy with one electron only. For two electrons with no repulsive Coulomb interaction, we have the following Schroedinger equation

~2 <em>d</em>2 <em>− </em>~2 <em>d</em>22 + 12<em>kr                     </em>2<em>kr </em>                                   <em>u</em>(<em>r</em>1<em>,r</em>2)<em>.</em>

<em>−             </em>2           2<em>m dr</em>2

2<em>m dr</em>1

Note that we deal with a two-electron wave function <em>u</em>(<em>r</em><sub>1</sub><em>,r</em><sub>2</sub>) and two-electron energy <em>E</em><sup>(2)</sup>.

With no interaction this can be written out as the product of two singleelectron wave functions, that is we have a solution on closed form.

We introduce the relative coordinate <strong>r </strong>= <strong>r</strong><sub>1 </sub><em>− </em><strong>r</strong><sub>2 </sub>and the center-of-mass coordinate <strong>R </strong>= 1<em>/</em>2(<strong>r</strong><sub>1</sub>+<strong>r</strong><sub>2</sub>). With these new coordinates, the radial Schroedinger equation reads

~2 <em>d</em>2 <em>− </em>~2 <em>d</em>2 + 1<em>kr</em>2 + <em>kR</em>2<em>u</em>(<em>r,R</em>) = <em>E</em>(2)<em>u</em>(<em>r,R</em>)<em>. − m dr</em><sup>2 </sup>4<em>m dR</em><sup>2 </sup>4

The equations for <em>r </em>and <em>R </em>can be separated via the ansatz for the wave function <em>u</em>(<em>r,R</em>) = <em>ψ</em>(<em>r</em>)<em>φ</em>(<em>R</em>) and the energy is given by the sum of the relative energy <em>E<sub>r </sub></em>and the center-of-mass energy <em>E<sub>R</sub></em>, that is

<em>E</em>(2) = <em>E</em><em>r </em>+ <em>E</em><em>R.</em>

We add then the repulsive Coulomb interaction between two electrons, namely a term

<em>V </em>(<em>r</em><sub>1</sub><em>,r</em><sub>2</sub>) = with <em>βe</em><sup>2 </sup>= 1<em>.</em>44 eVnm.




<em>βe</em><sup>2                      </sup><em>βe</em><sup>2</sup>

=          <em>,</em>

<em>|</em><strong>r</strong><sub>1 </sub><em>− </em><strong>r</strong><sub>2</sub><em>|                r</em>




Adding this term, the <em>r</em>-dependent Schroedinger equation becomes

~2 <em>d</em>2                  1     2             <em>βe</em>2

<em>− </em><sup>2 </sup>+ 4<em>kr </em>+ <em>r ψ</em>(<em>r</em>) = <em>E</em><em>rψ</em>(<em>r</em>)<em>. m dr</em>

This equation is similar to the one we had previously in (b) and we introduce again a dimensionless variable <em>ρ </em>= <em>r/α</em>. Repeating the same steps as above, we arrive at

<em>− </em><em>d</em><sup>2</sup>2 <em>ψ</em>(<em>ρ</em>) + 14 <em>mk</em>~2 <em>α</em><sub>4</sub><em>ρ</em><sub>2</sub><em>ψ</em>(<em>ρ</em>) + <em>mαβe</em><em>ρ</em>~2 <sup>2 </sup><em>ψ</em>(<em>ρ</em>) = <em>mα</em>~2<sup>2 </sup><em>E<sub>r</sub>ψ</em>(<em>ρ</em>)<em>. dρ</em>

We want to manipulate this equation further to make it as similar to that in (a) as possible. We define a new ’frequency’

2 1 <em>mk </em>4 <em>ω<sub>r </sub></em>= 4 ~<sub>2 </sub><em>α ,</em>

and fix the constant <em>α </em>by requiring

<em>mαβe</em><sup>2</sup>

= 1

~2

or ~2 <em>α </em>=                                                  <em>. mβe</em><sup>2</sup>

Defining

<em>mα</em><sup>2</sup>

<em>λ </em>=    <em>E,</em>

~2

we can rewrite Schroedinger’s equation as

<em>− </em><em>d</em>2<sup>2 </sup><em>ψ</em>(<em>ρ</em>) + <em>ω<sub>r</sub></em>2<em>ρ</em>2<em>ψ</em>(<em>ρ</em>) + <em>ρ</em>1 = <em>λψ</em>(<em>ρ</em>)<em>. dρ</em>

We treat <em>ω<sub>r </sub></em>as a parameter which reflects the strength of the oscillator potential.

Here we will study the cases <em>ω<sub>r </sub></em>= 0<em>.</em>01, <em>ω<sub>r </sub></em>= 0<em>.</em>5, <em>ω<sub>r </sub></em>= 1, and <em>ω<sub>r </sub></em>= 5 for the ground state only, that is the lowest-lying state.

With no repulsive Coulomb interaction you should get a result which corresponds to the relative energy of a non-interacting system. Make sure your results are stable as functions of <em>ρ</em><sub>max </sub>and the number of steps.

We are only interested in the ground state with <em>l </em>= 0. We omit the center-ofmass energy.

You can reuse the code you wrote for part (b), but you need to add the potential <em>ω<sub>r</sub></em><sup>2</sup><em>ρ</em><sup>2 </sup>+ 1<em>/ρ </em>to the diagonal elements.

Comment the results for the lowest state (ground state) as function of varying strengths of <em>ω<sub>r</sub></em>.

For specific oscillator frequencies, the above equation has answers in an analytical form, see the article by <a href="https://prola.aps.org/abstract/PRA/v48/i5/p3561_1">M. Taut, Phys. Rev. A 48, 3561 (1993)</a><a href="https://prola.aps.org/abstract/PRA/v48/i5/p3561_1">.</a>

<strong>Project 2 f): First optional exercise: Quantum physics analysis of the results. </strong>This exercise is a continuation of the previous and adds more quantum physics to the analysis. In this exercise we want to plot the wave function for two electrons as functions of the relative coordinate <em>r </em>and different values of <em>ω<sub>r</sub></em>. With no Coulomb interaction you should have a result which corresponds to the non-interacting case. Plot either the function itself or the probability distribution (the function squared) with and without the repulsion between the two electrons. Varying <em>ω<sub>r</sub></em>, the shape of the wave function will change.

We are only interested in the wave function for the ground state with <em>l </em>= 0 and omit again the center-of-mass motion.

The eigenvectors are normalized. Plot then the normalized wave functions for different values of <em>ω<sub>r </sub></em>and comment the results.

<strong>Project 2 g): Second optional challenge: Smarter way of finding the eigenvalues. </strong>This exercise is also optional and is an algorithmic challenge. Your matrix is already tridiagonal. Can you devise a more efficient way to find the eigenvalues and eigenvectors instead of the brute force application of Jacobi’s method? Hint: Use the bisection method discussed in for example <a href="http://www.maths.ed.ac.uk/~aar/papers/bamawi.pdf">Barth, Martin </a><a href="http://www.maths.ed.ac.uk/~aar/papers/bamawi.pdf">and Wilkinson, Numerische Mathematik 9, 386 (1967)</a><a href="http://www.maths.ed.ac.uk/~aar/papers/bamawi.pdf">.</a>

<strong>Project 2 h): Third optional challenge: Implementing Lanczos’ algorithm. </strong>This exercise is optional and is meant more as a challenge for those of you with an interest in other algorithms for solving eigenvalue problems. Implement the iterative Lanczos’ algorithm discussed in the lecture notes and compute the lowest eigenvalues as done in exercise (c) or (d) above. Compare your results and discuss faults and merits of the iterative method versus direct methods like Jacobi’s method.

<h1>Introduction to numerical projects</h1>

Here follows a brief recipe and recommendation on how to write a report for each project.

<ul>

 <li>Give a short description of the nature of the problem and the eventual numerical methods you have used.</li>

 <li>Describe the algorithm you have used and/or developed. Here you may find it convenient to use pseudocoding. In many cases you can describe the algorithm in the program itself.</li>

 <li>Include the source code of your program. Comment your program properly.</li>

 <li>If possible, try to find analytic solutions, or known limits in order to test your program when developing the code.</li>

 <li>Include your results either in figure form or in a table. Remember to label your results. All tables and figures should have relevant captions and labels on the axes.</li>

 <li>Try to evaluate the reliabilty and numerical stability/precision of your results. If possible, include a qualitative and/or quantitative discussion of the numerical stability, eventual loss of precision etc.</li>

 <li>Try to give an interpretation of you results in your answers to the problems.</li>

 <li>Critique: if possible include your comments and reflections about the exercise, whether you felt you learnt something, ideas for improvements and other thoughts you’ve made when solving the exercise. We wish to keep this course at the interactive level and your comments can help us improve it.</li>

 <li>Try to establish a practice where you log your work at the computerlab. You may find such a logbook very handy at later stages in your work, especially when you don’t properly remember what a previous test version of your program did. Here you could also record the time spent on solving the exercise, various algorithms you may have tested or other topics which you feel worthy of mentioning.</li>

</ul>

<h1>Format for electronic delivery of report and programs</h1>

The preferred format for the report is a PDF file. You can also use DOC or postscript formats or as an ipython notebook file. As programming language we prefer that you choose between C/C++, Fortran2008 or Python. The following prescription should be followed when preparing the report:

<ul>

 <li>Use Devilry to hand in your projects, log in at <a href="http://devilry.ifi.uio.no/">http://devilry.ifi. </a><a href="http://devilry.ifi.uio.no/">no</a> with your normal UiO username and password and choose either</li>

</ul>

’fys3150’ or ’fys4150’. There you can load up the files within the deadline.

<ul>

 <li>Upload <strong>only </strong>the report file! For the source code file(s) you have developed please provide us with your link to your github domain. The report file should include all of your discussions and a list of the codes you have developed. Do not include library files which are available at the course homepage, unless you have made specific changes to them.</li>

 <li>In your git repository, please include a folder which contains selected results. These can be in the form of output from your code for a selected set of runs and input parameters.</li>

 <li>In this and all later projects, you should include tests (for example unit tests) of your code(s).</li>

 <li>Comments from us on your projects, approval or not, corrections to be made etc can be found under your Devilry domain and are only visible to you and the teachers of the course.</li>

</ul>

Finally, we encourage you to work two and two together. Optimal working groups consist of 2-3 students. You can then hand in a common report.