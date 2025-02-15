(
var freq = 440;
var phase = 0;

// NB: iphase 0->1 is pi rads

{
    var sine, saw, saw2, triangle, square, square2;

    sine = SinOsc.ar(freq: freq, phase: phase * 2pi);
	//saw = Saw.ar(freq);
    saw = LFSaw.ar(freq: freq, iphase: phase);
	saw2 = LFSaw.ar(freq: freq, iphase: phase + 1);
    triangle = LFTri.ar(freq: freq, iphase: phase);
	square = Pulse.ar(freq: freq, width: 0.5);
    square2 = LFPulse.ar(freq: freq, iphase: phase + 1, width: 0.5);

	[sine, saw, saw2, triangle, square, square2]
}.plot(duration: 0.01); // 10ms
)
