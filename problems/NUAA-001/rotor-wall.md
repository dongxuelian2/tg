# NUAA-001 working manuscript — one-sided moving-wall viability with coupled rotor forces

Research draft, 2026-11-19. The result concerns finite-window safety in a planar rigid-body subsystem. It does not establish a theorem about arbitrary six-degree-of-freedom swarms.

## Dynamics and quantifiers

Two independently commanded rotor banks have normalized forces L,R in [0,1]. Set u=(R+L)/2, v=R-L. The exact mixer diamond is |v|/2 <= u <= 1-|v|/2. In a principal pitch plane with no yaw or roll torque the rigid-body equations are

    x'' = A u sin(theta) + w,   theta'' = B v,   z'' = A u cos(theta) - g.

A and B are positive constants. Initial x, horizontal velocity, pitch and pitch rate vanish. Altitude is unrestricted during this window (an open finite vertical band containing all possible trajectories works as well). Inputs are measurable; bank force slew is initially unrestricted. Disturbance w is measurable with |w|<=d. The wall constrains the center after accounting for a fixed body envelope. It spans the other spatial directions, so altitude alone cannot clear it. The planar restriction is part of the flight mode. We do not assert that allowing unrestricted three-dimensional attitude preserves the bound.

Fix H>0 with BH^2/2 <= pi/3. Let phi(t) denote the zero-wind position response to a specified pair of bank commands. Weak-star convergence of these commands in L-infinity implies uniform convergence of pitch, pitch rate, position and velocity. Indeed the pitch integrals converge pointwise and are equicontinuous; substitute uniform convergence of sin(theta) in the horizontal integrals and apply weak-star convergence to the remaining fixed kernels. The admissible controls form a sequentially weak-star compact set (bounded closed convex pointwise diamond, separable L1 predual). Thus finite-time maxima exist.

## Lemma — unique reset-free extremizer on an interval

On a subinterval of length h starting at global time b, suppose the incoming pitch and pitch rate are nonnegative and physically reachable from level rest. To maximize final x for zero wind, one can replace v by |v| and take u=1-|v|/2. The new pitch majorizes the absolute old pitch because the free pitch is nonnegative. All absolute pitches stay below pi/3 by the global bound, so the replacement pointwise majorizes horizontal acceleration. Incoming x and x' contribute an additive constant.

Write k=v in the reduced class, 0<=k<=1, and set U=1-k/2. For an optimum the sign of the first-variation density, divided by A(h-t), is the sign of

    D(t) = B integral_[t,h] [(h-s)(s-t)/(h-t)] U(s) cos(theta(s)) ds
           - sin(theta(t))/2.

This follows by differentiating along convex admissible variations and exchanging the two integrals. Needle variations, equivalently testing against arbitrary box-valued controls, give k=1 where D>0 and k=0 where D<0 almost everywhere. The integral decreases strictly in t: for s<h its kernel derivative is -(h-s)^2/(h-t)^2 and U cos(theta)>=1/4. The sine term is nondecreasing. Every optimum therefore has a prefix k=1 followed by k=0, including possible empty arcs. At the final time the integral tends to zero. For the first interval D(0)>0 and the terminal pitch of an optimum is positive. For later intervals with h>=3b we have

    D(0) >= B h^2/24 - B b^2/4 > 0

when b>0, using the incoming pitch <=Bb^2/2. In both cases the switch is strictly internal.

Uniqueness needs an additional argument since D depends on the trajectory. For a putative switch q, on the suffix s>=q the pitch is

    theta_q(s)=theta_in + omega_in s + B(q s-q^2/2).

The necessary scalar equation is

    B integral_[q,h] [(h-s)(s-q)/(h-q)] cos(theta_q(s)) ds
        = sin(theta_q(q))/2.

The left side decreases strictly with q by pointwise kernel comparison, shrinking integration domain, and monotonicity of theta_q(s) at fixed s. The right side increases strictly for q>0. Existence was proved by compactness and the first variation, so the root is unique, and bisection constructs it. In the unreduced controls equality in the majorization is necessary almost everywhere. The reduced optimum has U>=1/2 and strictly positive pitch for t>0. Equality forces the same u and the same pitch a.e. (positive sine and strict monotonicity), hence the same v by twice differentiating the equal pitch paths. This proves uniqueness up to null sets.

## Theorem — arbitrary control variation forced by one moving wall

For every integer n>=1, fixed H as above, and any d>=0, there exists one C^(1,1) moving half-space x>=f(t), t in [0,H], with a robustly feasible positive clearance. Nevertheless every robustly safe controller on its constant headwind realization has normalized differential force of total variation greater than 2n-2 (essential variation, independent of null-set modifications). Wall speed and acceleration bounds can be chosen independently of n. In particular no bound on the number of constant-command segments follows from horizon, vehicle force limits, and these wall derivative bounds alone. The clearance is allowed to depend on n.

Proof. Choose t_j=H 4^(j-n), 1<=j<=n, t_0=0. Successively maximize x(t_j) from the preceding maximizer's terminal state, with zero wind. After the first interval each duration equals three times its start time, so the lemma gives an internal switch on every interval. Let X denote the concatenated position. In every segment v alternates 1,0. The concatenation has essential variation 2n-1. Among zero-wind controls, satisfying x(t_j)>=X(t_j) for every j forces exactly this concatenation. Apply uniqueness first at j=1 and then inductively to each fixed incoming state. For nonzero wind the same argument applies on the constant realization w=-d after adding dt^2/2 to the realized position.

For a strict margin define the finite scalar gap

    G_n = min_[TV(v)<=2n-2] max_j (X(t_j)-phi(t_j)).

The minimum exists: bounded-variation v has an L1-convergent subsequence by Helly compactness, u has a weak-star convergent subsequence, the diamond is closed, and variation is lower semicontinuous. Every member has a strictly positive maximum deficit by uniqueness, hence G_n>0. Use the wall f(t)=X(t)-dt^2/2-epsilon with 0<epsilon<G_n. The concatenated open-loop controls guarantee clearance epsilon continuously, since the worst horizontal displacement due to wind is -dt^2/2. Any purported robust controller with headwind TV<=2n-2 violates at least one checkpoint. Feedback cannot alter the exact pathwise input-output relation. The wall has |f'|<=(A+d)H and |f''|<=A+d a.e., independently of n. With d=0 it is also nondecreasing and convex: the obstacle need never oscillate or reverse. A constant-command signal with k switches has TV<=2k since |v|<=1; thus k>=n is necessary. More directly compactness of the family with <=2n-2 switches, with free command amplitudes, and the unique 2n-1-switch tight trajectory give the sharper lower bound 2n-1 switches for sufficiently small epsilon. This last bound follows by using the minimum of the two positive gaps. QED.

The strict margin also allows smoothing the wall by uniform approximation to error less than half of both margins (feasibility and exclusion); linear tangent extensions outside [0,H] followed by convolution preserve uniform first and second derivative bounds. For each fixed n the witnessing rotor steps can be smoothed by short ramps inside the diamond with uniformly small position error. Thus the conclusion survives a sufficiently large finite slew limit for each instance, with a feasible rate-limited witness. We have not shown an n-independent slew bound.

## An information obstruction with positive margin

Consider two environments with mirrored walls and safe sides, using the same parameters and reflected witness (v changes sign). A remote observer knows the environment; the executor has identical initial information and receives no distinguishing information until a prescribed delta in (0,H]. The common-prefix assumption includes its onboard observations on the zero-wind realizations. This can model a hidden gate signal; it would not describe an executor that can see the wall immediately. No communication is needed after a single correct initial bit.

At d=epsilon=0 the two feasible paths are individually unique, and their v have opposite signs on a common initial interval of positive length. Pairs of admissible zero-wind controls equal a.e. on [0,delta] form a weak-star compact set. Minimize over this set the maximum of all signed checkpoint deficits in both worlds. Call the value D_delta. By the same uniqueness argument D_delta>0. For any epsilon,d>0 with epsilon+dH^2/2<D_delta, at least one zero-wind execution violates a checkpoint even if perfect information is delivered at delta. Immediate one-bit disclosure is sufficient under every allowed wind. Taking epsilon<G_n simultaneously gives the variation obstruction. This proves a strict latency separation for arbitrary causal flight controllers under the stipulated information pattern. It gives no general bandwidth law when observations of the wall are available earlier.

## Relatedness and open audit

* [Guo et al., flexible virtual tubes](https://doi.org/10.1016/j.cja.2024.11.010) use positive tube certificates. This manuscript addresses exclusion of entire control classes in a single projected rigid-body safety task.
* [Quadcopter position isochrones](https://par.nsf.gov/servlets/purl/10318015) use a bounded angular velocity input. An initial angular acceleration and its consumption of collective thrust are essential here.
* [Foehn et al.](https://rpg.ifi.uzh.ch/docs/ScienceRobotics21_Foehn.pdf) explicitly account for dependent thrust/torque resources. The existence of the mixer coupling itself is prior knowledge.
* [Qin et al., IJRR 2026](https://doi.org/10.1177/02783649261480153) analyze point-to-point time optimality and discuss extending to intermediate waypoints. The moving-wall theorem is not a counterexample to their point-to-point switch bounds.
* [Caponigro et al.](https://arxiv.org/abs/1303.5796) study variation regularization and chattering. The theorem here gives finite-window robust feasibility separation for one monotone wall in its disturbance-free specialization; no single instance of infinite chattering is claimed.

Open issues include quantitative estimates for the robustness gaps as n increases, a genuinely coupled multi-vehicle obstruction, matching feasibility conditions with finite motor bandwidth fixed across the family, and an expanded priority search in optimal-control and viability literature. A repository path for the user's TiGuan record has not yet been verified. No publication priority or completion claim.
