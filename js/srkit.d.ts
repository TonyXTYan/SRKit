// Type declarations for srkit.js: SRKit for JavaScript.
// Natural units, c = 1: one unit of time, and the distance light covers in it as the unit of length.

export type Vec3 = [number, number, number];

/** At rest at P. */
export interface RestSegment { kind: 'rest'; t0: number; t1: number; P: Vec3; tau0?: number }
/** Constant velocity V (a fraction of c) from P0. `ig` is 1 / gamma. */
export interface LinSegment { kind: 'lin'; t0: number; t1: number; P0: Vec3; V: Vec3; ig?: number; tau0?: number }
/** Proper acceleration `al` (signed) along the unit vector n; the rapidity along n is ph0 at t0. */
export interface HypSegment { kind: 'hyp'; t0: number; t1: number; P0: Vec3; n: Vec3; al: number; ph0: number; tau0?: number }
/** A circle about C of radius r in the plane (e1, e2), at angle a0 + w (t - t0). */
export interface CircSegment { kind: 'circ'; t0: number; t1: number; C: Vec3; r: number; w: number; e1: Vec3; e2: Vec3; a0: number; ig?: number; tau0?: number }
/** A light-like jump from P0 along n to P1. No proper time passes. */
export interface NullSegment { kind: 'null'; t0: number; t1: number; P0: Vec3; n: Vec3; P1: Vec3; tau0?: number }
/** Proper acceleration `al` > 0 at any angle: `thrust` is a direction fixed aboard, U0 the starting proper velocity (gamma v). */
export interface SteerSegment { kind: 'steer'; t0: number; t1: number; P0: Vec3; U0: Vec3; thrust: Vec3; al: number; tau0?: number }
/** Constant thrust as the mass falls: starting rapidity ph0 along n, starting proper acceleration al0 (signed), exhaust speed w. */
export interface BurnSegment { kind: 'burn'; t0: number; t1: number; P0: Vec3; n: Vec3; ph0: number; al0: number; w: number; tau0?: number }
export type Segment = RestSegment | LinSegment | HypSegment | CircSegment | NullSegment | SteerSegment | BurnSegment;

export interface State { X: Vec3; V: Vec3 }
export interface TauState { t: number; X: Vec3; V: Vec3; phi: number; acc: number; jump?: true }

export interface SegmentApi {
  mk: {
    rest(t0: number, t1: number, P: Vec3): RestSegment;
    lin(t0: number, t1: number, P0: Vec3, V: Vec3): LinSegment;
    hyp(t0: number, t1: number, P0: Vec3, n: Vec3, al: number, ph0: number): HypSegment;
    circ(t0: number, t1: number, C: Vec3, r: number, w: number, e1: Vec3, e2: Vec3, a0: number): CircSegment;
    null(t0: number, t1: number, P0: Vec3, n: Vec3): NullSegment;
    steer(t0: number, t1: number, P0: Vec3, U0: Vec3, thrust: Vec3, al: number): SteerSegment;
    burn(t0: number, t1: number, P0: Vec3, n: Vec3, ph0: number, al0: number, w: number): BurnSegment;
  };
  /** Position and velocity at coordinate time t. */
  at(g: Segment, t: number): State;
  /** Proper time from t0 to t. */
  dtau(g: Segment, t: number): number;
  /** Coordinate time after proper time s from t0. */
  tOf(g: Segment, s: number): number;
  /** The state after proper time s from t0. */
  atTau(g: Segment, s: number): TauState;
  /** The largest speed on the segment. */
  vmax(g: Segment): number;
  /** When light from the event (tE, XE) reaches the segment, within [lo, hi]; null if it does not. */
  recv(g: Segment, tE: number, XE: Vec3, lo: number, hi: number): number | null;
  /** The time on the segment, within [lo, hi], whose light reaches the event; null if there is none. */
  send(g: Segment, tE: number, XE: Vec3, lo: number, hi: number): number | null;
  /** For a hyp or steer segment continued for ever: how far the event lies inside its horizon (positive: inside). */
  horizonMargin(g: Segment, tE: number, XE: Vec3, past?: boolean): number | null;
  properDuration(g: Segment): number;
}

/** A worldline as a plain object. `make` sets tau0 on every segment and has the engine hold a copy. */
export interface WorldlineObject { segs: Segment[]; t0s?: number[]; tau0s?: number[]; vmax?: number; hasNull?: boolean }
export interface WorldlineApi {
  make(segs: Segment[], anchorT?: number, anchorTau?: number): WorldlineObject;
  find(W: WorldlineObject, t: number): number;
  findTau(W: WorldlineObject, tau: number): number;
  at(W: WorldlineObject, t: number): State;
  tau(W: WorldlineObject, t: number): number;
  atTau(W: WorldlineObject, tau: number): TauState;
  /** When light from the event reaches the worldline; Infinity if never. */
  recv(W: WorldlineObject, tE: number, XE: Vec3): number;
  /** The time on the worldline whose light reaches the event; -Infinity if none. */
  send(W: WorldlineObject, tE: number, XE: Vec3): number;
  options: { handles: boolean };
  /** Worldlines the engine currently holds. */
  live(): number;
}

/** A worldline owned by the engine. Release it with free(), or leave it to the garbage collector. */
export declare class Worldline {
  static from(segs: Segment[], anchorT?: number, anchorTau?: number): Worldline;
  static fromAnchored(segs: Segment[]): Worldline;
  readonly length: number;
  free(): void;
  at(t: number): State & { tau: number };
  tau(t: number): number;
  tOfTau(tau: number): number;
  atTau(tau: number): TauState;
  find(t: number): number;
  recv(tE: number, XE: Vec3): number;
  send(tE: number, XE: Vec3): number;
  vmax(): number;
  /** Whether consecutive segments join in position and velocity. */
  continuous(positionTol?: number, velocityTol?: number): boolean;
  /** Radar coordinates of the event: the proper times at which a ping must leave and its echo returns. */
  radar(tE: number, XE: Vec3): { tauSend: number; tauRecv: number };
  segment(i: number): Segment | null;
  /** Seven numbers per time: x, y, z, vx, vy, vz, tau. */
  atBatch(ts: ArrayLike<number>): Float64Array;
  /** `events` is flat: t, x, y, z per event. Returns the time on this worldline for each. */
  coneBatch(events: ArrayLike<number>, send?: boolean): Float64Array;
}

export interface Start { t: number; tau: number; X: Vec3 }
export interface End { t: number; tau: number; X: Vec3; V: Vec3; phi: number; acc: number }
export type Phase = 'rest' | 'boost' | 'coast' | 'brake' | 'orbit' | 'jump';
export interface Piece { g: Segment; dtau: number; phase: Phase }
export interface Line { pieces: Piece[]; end: End; n: Vec3; dphi: number; L: number; phPeak: number; coastL: number }
export interface OrbitElements { C: Vec3; r: number; w: number; e1: Vec3; e2: Vec3; a0: number; beta?: number }
export interface PlanApi {
  rest(start: Start, dtau: number): Piece;
  hypPiece(start: Start, n: Vec3, al: number, ph0: number, ph1: number, phase: Phase): Piece;
  linPiece(start: Start, V: Vec3, dtau: number): Piece;
  nullPiece(start: Start, X1: Vec3): Piece;
  circPiece(start: Start, orb: OrbitElements, dtau: number): Piece;
  pieceEnd(p: Piece): TauState & { tau: number };
  /** A straight run: boost, coast if capped at rapidity phCap, brake. null if too short for the speeds asked. */
  line(start: Start, X1: Vec3, a: number, phCap?: number, ph0?: number, ph1?: number, aBrake?: number, nDefault?: Vec3 | null): Line | null;
  brake(start: Start, n: Vec3, ph0: number, a: number): { pieces: Piece[]; end: End; dphi: number };
  orbitExit(start: Start, V: Vec3, a: number): { pieces: Piece[]; end: End; dphi: number };
  planeThrough(C: Vec3, X: Vec3, e1Default?: Vec3 | null): { e1: Vec3; e2: Vec3; nrm: Vec3 };
  orbitEntry(start: Start, C: Vec3, r: number, beta: number, a: number, phCap?: number, opts?: { e1?: Vec3; a0?: number; aTangent?: number }):
    { pieces: Piece[]; end: End; orb: OrbitElements; dphi: number; fromInside: boolean } | null;
}

/** A ride of rest-to-rest legs in null coordinates (u = t - x, v = t + x). */
export interface Ride {
  d: number; a: number; n: number; nOut: number; rt: boolean; legTau: number; tau: number; du: number; dv: number;
  uEnd: number; vEnd: number; tEnd: number; tOut: number; xEnd: number; lag: number; gPeak: number; phc: number;
  leg: { a: number; phc: number; ta: number; tc: number };
}
export interface RideState { tau: number; u: number; v: number; t: number; x: number; phi: number; dir: number; gamma: number; Dh: number; Dd: number }
export interface TripApi {
  plan(d: number, a: number, cap?: number, stops?: number, rt?: boolean): Ride;
  /** The ship at proper time tau; it waits at rest before launch and after arrival. */
  at(p: Ride, tau: number): RideState;
  /** As `at`, held within the ride. */
  atClamped(p: Ride, tau: number): RideState;
  tauU(p: Ride, u: number): number;
  tauV(p: Ride, v: number): number;
  /** Proper time at which the ride meets a light cone of the event (t, along, across). */
  cone(p: Ride, t: number, along: number, across: number, send?: boolean): number;
  flipBurn(d: number, a: number): { d: number; a: number; t: number; tau: number; phiPeak: number; gammaPeak: number; lag: number };
  /** Where a ship at rapidity phi can stop within a proper-time budget, burning at `heading` from its velocity. */
  reach(phi: number, a: number, budget: number, heading: number): [number, number] | null;
}

export interface Orbit { r: number; beta: number; gamma: number; phi: number; omega: number; period: number; lapTau: number; acc: number; sagnac: number; thomas: number }
export interface OrbitApi {
  ofSpeed(r: number, beta: number): Orbit;
  ofAcceleration(r: number, acc: number): Orbit;
  at(r: number, acc: number, a0: number, tau: number): { t: number; X: Vec3; V: Vec3; heading: Vec3; al: number };
  cone(r: number, acc: number, a0: number, t: number, X: Vec3, send?: boolean): number;
}

export interface OpticsApi {
  /** Aberration and Doppler for a source in direction (a, br, bu) in the observer's forward, right and up axes. */
  observe(phi: number, a: number, br: number, bu: number, aberration?: boolean): { D: number; th: number; thp: number; f: number; r: number; u: number };
  playback(V: Vec3, n: Vec3, VB: Vec3): number;
  basis(V: Vec3): { f: Vec3; u: Vec3; r: Vec3 };
  /** Boost (dt, dX) into the frame whose velocity points along w with rapidity |w|. */
  boost(w: Vec3, dt: number, dX: Vec3): { t: number; X: Vec3 };
  compose(v: Vec3, u: Vec3): Vec3;
  resolution(aperture: number, wavelength: number, distance: number): number;
  friis(power: number, dT: number, dR: number, wavelength: number, distance: number): number;
  bitRate(powerW: number, wavelengthM: number, bitsPerPhoton?: number): number;
  luminosity(radiusM: number, temperatureK: number): number;
}

export interface RocketApi {
  massRatio(dphi: number, w?: number): number;
  phiOf(massRatio: number, w?: number): number;
  /** Propellant energy per unit of dry rest energy. Multiply by m c² for joules. */
  energy(dphi: number, w?: number, efficiency?: number): number;
  sailPhi(energyPerRestEnergy: number): number;
}

/** Whole milliseconds as decimal strings or BigInt; results come back in the type of the first argument. */
export type Exact = string | bigint;
export interface ExactTime {
  bits: 384 | 512;
  add<T extends Exact>(t: T, d: Exact): T;
  diff<T extends Exact>(a: T, b: Exact): T;
  cmp(a: Exact, b: Exact): -1 | 0 | 1;
  /** The exact phase of t in a period of whole milliseconds. */
  mod(t: Exact, periodMs: number): number | null;
  /** a - b in seconds, subtracted exactly first; null beyond 2^53 ms. */
  seconds(a: Exact, b: Exact): number | null;
  /** Light-travel time over dX metres, rounded up to whole milliseconds. */
  lightTicks(dX: Vec3): string;
  firstArrival<T extends Exact>(emission: T, dX: Vec3): T;
  /** Has a signal emitted at `emission` crossed dX metres by the tick `now`? Exact. */
  arrived(emission: Exact, now: Exact, dX: Vec3): boolean;
  classify(dt: Exact, dX: Vec3): 'timelike' | 'null' | 'spacelike';
}

/** A flat universe of matter and a cosmological constant, seen from inside. Times in years since the Big Bang, lengths in
 *  light-years; comoving lengths are proper lengths today. Radiation is not included. */
export interface Cosmology {
  h0: number; omegaM: number;
  age(): number; scaleFactor(t: number): number; timeAtScaleFactor(a: number): number; hubble(t: number): number;
  conformalTime(t: number): number; conformalTimeLimit(): number;
  eventHorizon(t: number): number; eventHorizonProper(t: number): number; particleHorizon(t: number): number; particleHorizonProper(t: number): number;
  /** How far (comoving) light sent at t can ever get. */
  reach(t: number): number;
  /** The last time a signal can be sent to a target at comoving distance chi and still arrive; null if it is already too late. */
  lastContactTime(chi: number): number | null;
  /** The last time to send and still receive a reply. */
  lastReplyTime(chi: number): number | null;
  /** Light cones in an expanding universe; null if the light never arrives. */
  lightArrivalTime(tEmitted: number, chi: number): number | null;
  lightEmissionTime(tReceived: number, chi: number): number | null;
  comovingDistance(tFrom: number, tTo: number): number; lookbackTime(z: number): number; timeAtRedshift(z: number): number; redshiftBetween(tEmitted: number, tReceived: number): number;
  /** Hubble drag: gamma beta of a coasting body falls as 1 / a. */
  coastingProperVelocity(u: number, tFrom: number, tTo: number): number;
  /** The Lorentz factor at which a rocket's thrust (alpha per year) balances Hubble drag. */
  dragLimitedGamma(alpha: number): number;
  hLambda(): number; deSitterHorizon(): number;
  /** Kelvin. */
  deSitterTemperature(): number; cmbTemperature(t: number, tCmb0?: number): number; cmbCrossoverTime(tCmb0?: number): number;
  /** Kilograms within reach of light sent at t. */
  reachMass(t: number): number; reachableShare(t: number): number;
}

/** SI throughout: kg, m, s, K, W. A spin is the dimensionless a*; where a function takes only a spin, negative means a retrograde orbit. */
export interface GravityApi {
  /** In units of GM/c². */
  kerrIscoRadius(spin: number): number; kerrHorizonRadius(spin: number): number;
  kerrIscoDtDtau(spin: number): number; accretionEfficiency(spin: number): number; extractableEnergyFraction(spin: number): number;
  zamoDtauDt(massKg: number, spin: number, r: number, polarAngle?: number): number;
  schwarzschildRadius(massKg: number): number; staticDtauDt(massKg: number, r: number): number; circularOrbitDtDtau(massKg: number, r: number): number;
  redshift(massKg: number, rEmitter: number, rObserver: number): number;
  shapiroDelay(massKg: number, rEmitter: number, rReceiver: number, impactParameter: number): number;
  deflectionAngle(massKg: number, impactParameter: number): number;
  hawkingTemperature(massKg: number): number; hawkingPower(massKg: number): number; hawkingEntropy(massKg: number): number;
  /** 'photons': a black body over the horizon area (the oracle's model). 'page': Page 1976, photons and gravitons. */
  hawkingLifetime(massKg: number, model?: 'photons' | 'page'): number; hawkingMassForLifetime(seconds: number): number;
}

export interface Units {
  C: number; YEAR: number; LY: number; AU: number; PC: number; G0: number; H: number; HBAR: number; KB: number; G: number; SIGMA: number;
  L_SUN: number; R_SUN: number; GM_SUN: number; G_PER_YEAR: number; AU_LY: number; PC_LY: number;
}

export interface SRKitApi {
  info: { backend: 'srkit-wasm'; version: string; abi: number; sha256: string; bytes: number };
  /** Calls made into the WebAssembly module so far. */
  calls(): number;
  units: Units;
  seg: SegmentApi;
  wl: WorldlineApi;
  Worldline: typeof Worldline;
  plan: PlanApi;
  trip: TripApi;
  orbit: OrbitApi;
  optics: OpticsApi;
  rocket: RocketApi;
  kepler(meanAnomaly: number, e: number): number;
  time: ExactTime & { u384: ExactTime; u512: ExactTime };
  /** The Planck 2018 universe (srkit-flrw); `flrw.cosmology(h0PerYear, omegaM)` makes another. */
  flrw: Cosmology & { cosmology(h0PerYear?: number, omegaM?: number): Cosmology; fromKmSMpc(h0KmSMpc: number, omegaM: number): Cosmology };
  /** Closed-form general-relativity scalars (srkit-gr). */
  gr: GravityApi;
  /** The module's raw exports. */
  raw: Record<string, unknown>;
  /** Instantiate another copy from WebAssembly bytes. */
  create(bytes: BufferSource): SRKitApi;
}

declare const SRKit: SRKitApi;
export default SRKit;
export as namespace SRKit;
