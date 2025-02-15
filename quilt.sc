
// issues:
// - TZO (ring mod) not working as expected
// - `counter*` waves are not 0-aligned!
// - phase shift impl (sync w/ wip code on norns)

(
s.waitForBoot {

var win, freqSlider, ampSlider, freqLabel, ampLabel,
  mod1Label, mod1Slider, mod1AmtLabel, mod1AmtSlider, mod1FreqLabel, mod1FreqSlider,
  mod2Label, mod2Slider, mod2AmtLabel, mod2AmtSlider, mod2FreqLabel, mod2FreqSlider,
  i1Label, i1Slider, i2Label, i2Slider, i3Label, i3Slider, i4Label, i4Slider;
var ly, lh = 30;
var def;

var d_amp_offset = 1;
var d_freq = 200;
var d_mod1 = 2;
var d_mod2 = 1;
var d_mod1a = 0;
var d_mod2a = 0;
var d_mod1f = 1;
var d_mod2f = 20;

var g_freq = d_freq;
var g_mod1 = d_mod1;


~sawBuffer = Buffer.alloc(s, 4096, 1);
~sawValues = (0..(~sawBuffer.numFrames-1)).collect { |i| (i / (~sawBuffer.numFrames-1)) * 2 - 1 };
~sawBuffer.loadCollection(~sawValues);

// ------------------------------------
// helper fns

~updateScopeRange = { |freq, periods|
	var cycles = periods * (s.sampleRate / freq).asInteger;
	var adjustedFreq = s.sampleRate/(cycles/periods) ;

	Stethoscope.ugenScopes.do({
		arg scope, i;
		scope.cycle = cycles;
	});

	adjustedFreq;
};

~hzToVolts = { |freq|
	(freq / 20).log2;
};

~voltsToHz = { |volts|
	20 * (2 ** volts);
};

~instantCutoff = { |baseCutoffHz, cutoffOffnessPct, keyHz, keyTrackPct, keyTrackNegOffsetPct, eg, envelopePct|
			var baseCutoffVolts = ~hzToVolts.(baseCutoffHz);
			var cutoffOffnessVolts = cutoffOffnessPct * 10 / 4; // NB: max offness by 2.5 volt

			var keyVolts = ~hzToVolts.(keyHz);
			var keyTrackNegOffsetVolts = keyTrackNegOffsetPct * 10;
			var egVolts = eg * 10;

			var keyModVolts = (keyVolts - keyTrackNegOffsetVolts) * keyTrackPct;
			var egModVolts  = egVolts  * envelopePct;

			var totalVolts = baseCutoffVolts + cutoffOffnessVolts + keyModVolts + egModVolts;

			var instantCutoffHz = ~voltsToHz.(totalVolts);

			instantCutoffHz = instantCutoffHz.clip(20, 20000);

			instantCutoffHz;
};


// ------------------------------------

		def = SynthDef(\Quilt, {
			arg out = 0,
			gate = 0,
			gate_pair_in = 0,
			gate_pair_out = 0,
			vel = 0.5,
			freq = 200,
			raw_osc_cutoff = 10000,
			phased_cutoff = 10000,
			freq_sag = 0.1,
			vib_rate = 5,
			vib_depth = 0.0,
			// segmented oscilators
			index1 = 0.0,
			index2 = 0.0,
			index3 = 0.0,
			index4 = 0.0,
			amp1 = 0.5,
			amp2 = 0.5,
			amp3 = 0.5,
			amp4 = 0.5,
			// npolar projection
			mod = 3,
			syncRatio = 1,
			syncPhase = 0.0,
		    syncPhaseSliced = 0.0,
			npolarProj = 1.0,
			npolarRotFreq = 10,
			npolarRotFreq_sag = 0.1,
			npolarProjSliced = 1.0,
			npolarRotFreqSliced = 10,
			npolarRotFreqSliced_sag = 0.1,
			// phase mod
			pmFreq = 0.5,
			pmAmt = 0,
			// amp env
			amp_offset = 0.0,
			attack = 0.1, decay = 0.1, sustain = 0.7, release = 0.5,
			// filter env
			fenv_a = 1.0,
			fktrack = 0.1,
			fktrack_neg_offset = 0.0,
			fattack = 0.1, fdecay = 0.1, fsustain = 0.7, frelease = 0.5,
			// filter
			cutoff = 1200,
			cutoff_sag = 0.1,
			resonance = 0.0,
			// panning
			pan = 0,
			pan_lfo_amount = 0.1,
			pan_lfo_freq = 5,
			pan_lfo_phase = 0,
			// offness
			phase_offset = 0.0,
			pitch_offness_max = 0.0,
			pitch_offness_pct = 0.0,
			cutoff_offness_max = 0.0,
			cutoff_offness_pct = 0.0,
			// saturation/compression
			sat_threshold = 0.5
		;

			// frequencies
			var semitoneDiff, hzTrack, fsemitoneDiff;
			var vibrato, freqSagLfo, freq2;
			var cutoffSagLfo, cutoff2;
			var npolarRotFreqSagLfo, npolarRotFreq2;
			var npolarRotFreqSlicedSagLfo, npolarRotFreqSliced2;
			// basic waveforms
			var sin, saw, triangle, square;
			// looked-up waveforms (by index)
			var signal1, signal2, signal3, signal4;
			// amp enveloppe
			var env, pairingInEnv, pairingOutEnv, scaledEnv;
			// filter enveloppe
			var fenv, instantCutoff;
			// processed waveform
			var mixed, phased, filtered, ironed, saturated, compressed;

			// CMOS-derived waveforms
			var crossing, counter, crossingSliced, counterSliced;
			// computed modulation index, associated phaser signals
			var phaseAm, phaseRm, phaseRmFade, phaseAmFade, phase;
			var amToRm, amToRmSliced;
			var phaseSlicedAm, phaseSlicedRm, phaseSlicedRmFade, phaseSlicedAmFade, phaseSliced;
			var pm;

			vibrato = SinOsc.kr(vib_rate, 0, vib_depth);

			// NB: this sounds meh and is heavy in processing...
			// TODO: implement standard vibrrato w/ slightly detuned voices

			semitoneDiff = freq * (2 ** (1/12) - 1);
			fsemitoneDiff = cutoff * (2 ** (1/12) - 1);
			// freqSagLfo = Lag.kr(LFNoise1.kr(1), 0.1) * freq_sag * semitoneDiff;
			// freq2 = freq + freqSagLfo + vibrato;

			// cutoffSagLfo = Lag.kr(LFNoise1.kr(1), 0.1) * cutoff_sag * semitoneDiff;
			// cutoff2 = cutoff + cutoffSagLfo;

			// npolarRotFreqSagLfo = Lag.kr(LFNoise1.kr(1), 0.1) * npolarRotFreq_sag * semitoneDiff;
			// npolarRotFreq2 = npolarRotFreq + npolarRotFreqSagLfo;

			// npolarRotFreqSlicedSagLfo = Lag.kr(LFNoise1.kr(1), 0.1) * npolarRotFreqSliced_sag * semitoneDiff;
			// npolarRotFreqSliced2 = npolarRotFreqSliced + npolarRotFreqSlicedSagLfo;

			freq2 = freq + vibrato + ((semitoneDiff) * pitch_offness_max * pitch_offness_pct);
			// cutoff2 = cutoff + (7000 * cutoff_offness_max * cutoff_offness_pct);
			npolarRotFreq2 = npolarRotFreq;
			npolarRotFreqSliced2 = npolarRotFreqSliced;

			hzTrack = freq2.cpsmidi / 12;

			sin = SinOsc.ar(freq2) * 0.5; // FIX: needed to half amp for sine
			saw = MoogFF.ar(in: Saw.ar(freq2), freq: raw_osc_cutoff);
			triangle = MoogFF.ar(in: LFTri.ar(freq2), freq: raw_osc_cutoff);
			square = MoogFF.ar(in: Pulse.ar(freq: freq2, width: 0.5), freq: raw_osc_cutoff);

			pm = SinOsc.ar(pmFreq) * pmAmt;

			crossing = Osc.ar(~sawBuffer, freq2 * 2, pi + (pm + syncPhase).linlin(-1, 1, -2pi, 2pi)) * 0.25;
			counter = PulseCount.ar(crossing) % mod;

			crossingSliced = Osc.ar(~sawBuffer, freq2 * syncRatio * 2, pi + (pm + syncPhaseSliced).linlin(-1, 1, -2pi, 2pi)) * 0.25;
			counterSliced = PulseCount.ar(crossingSliced) % mod;

			// REVIEW: use wavetable instead?
			signal1 = Select.ar(index1, [sin, triangle, saw, square]);// * amp1 * SinOsc.kr(npolarRotFreq, 0.0);
			signal2 = Select.ar(index2, [sin, triangle, saw, square]);// * amp2 * SinOsc.kr(npolarRotFreq, 2pi / mod);
			signal3 = Select.ar(index3, [sin, triangle, saw, square]);// * amp3 * SinOsc.kr(npolarRotFreq, 2 * 2pi / mod);
			signal4 = Select.ar(index4, [sin, triangle, saw, square]);// * amp4 * SinOsc.kr(npolarRotFreq, 3 * 2pi / mod);

			mixed = Select.ar(counterSliced, [signal1, signal2, signal3, signal4]) * 2;

			phaseRm = SinOsc.ar(npolarRotFreq2, counter * 2pi/mod, 1);
		    //phaseAm = if(mod % 2 == 0, { phaseRm }, { (1.0 - phaseRm) }) / 2;
			// NB: edge-case for when mod1 is 2
			// this works, but idk why using `if(mod == 2, ...)` doesn't
			phaseRm = Select.ar((mod-2).clip(0, 1),
				[
			//SinOsc.ar(npolarRotFreq2, counter * 2pi/(mod-1), 1),
			-1 * phaseRm,
					phaseRm ]);
			// NB: phaseAmFade crossfades between a DC of 1 and phaseAm according to npolarProj
			// the critical part is the 3rd argument of phaseRmFade
			// i don't fully understand how it works...
		    phaseRmFade = SinOsc.ar(npolarRotFreq2, counter * 2pi/mod, (npolarProj*2).clip(0,1));
			phaseAmFade = if(mod % 2 == 0, { phaseRmFade }, { (1.0 - phaseRmFade) });

			// x-fade between phaseAmFade and phaseRm, according to npolarProj (only for 0.5-1)
			// we could have used XFade2.ar instead...
		    amToRm = (npolarProj-0.5).clip(0, 0.5) * 2;
			phase = (phaseAmFade * (1 - amToRm))
			+ (phaseRm * 2 * amToRm * (-1))
			;

			phaseSlicedRm = SinOsc.ar(npolarRotFreqSliced2, counterSliced * 2pi/mod, 1);
			//phaseSlicedAm = if(mod % 2 == 0, { phaseSlicedRm }, { (1.0 - phaseSlicedRm) }) / 2;
			phaseSlicedRmFade = SinOsc.ar(npolarRotFreqSliced2, counterSliced * 2pi/mod, npolarProjSliced);
			phaseSlicedAmFade = if(mod % 2 == 0, { phaseSlicedRmFade }, { (1.0 - phaseSlicedRmFade) });

			amToRmSliced = (npolarProjSliced-0.5).clip(0, 0.5) * 2;
			phaseSliced = (phaseSlicedAmFade * (1 - amToRmSliced))
			+ (phaseSlicedRm * 2 * amToRmSliced * (-1))
			;

			// phased = mixed
			// * ((npolarProj*2).linlin(0, 1, 1, phase))
			// * ((npolarProjSliced*2).linlin(0, 1, 1, phaseSliced));

			phased = mixed * phase * phaseSliced;

			phased =  MoogFF.ar(in: phased, freq: phased_cutoff, gain: 0);

			env = EnvGen.kr(Env.adsr(attack, decay, sustain, release), gate, doneAction: 0);
			// NB: enveloppes for when a voice is dynamically paired
			pairingInEnv  = EnvGen.kr(Env.adsr(0.7, 0, sustain, 0), gate_pair_in, doneAction: 0);
			pairingOutEnv = EnvGen.kr(Env.adsr(0, 0, sustain, 0.7), gate_pair_out, doneAction: 0);
			scaledEnv = (1 - amp_offset) * (env + pairingInEnv + pairingOutEnv) + amp_offset;

			// fenv = EnvGen.kr(Env.adsr(fattack, fdecay, fsustain, frelease), gate, doneAction: 0) * (fenv_a / 2);
			fenv = EnvGen.kr(Env.adsr(fattack, fdecay, fsustain, frelease), gate, doneAction: 0);

			// instantCutoff = (cutoff2 + (fktrack * (freq2.cpsmidi).clip(21, 127).linexp(21, 127, 27.5, 12543.85)) + fenv.linlin(0, 1, 0, 15000)).clip(20, 20000);
			instantCutoff = ~instantCutoff.(cutoff, cutoff_offness_max * cutoff_offness_pct,
				freq2, fktrack, fktrack_neg_offset,
				fenv, fenv_a);

			filtered = MoogFF.ar(in: phased,
				freq: instantCutoff,
				gain: resonance) * 0.5 * vel * scaledEnv;

			ironed = BPeakEQ.ar(filtered, 200, rq: 1, db: 6 * (1-sat_threshold));

			saturated = (ironed * (2 - sat_threshold)).tanh;

			// compressed = Compander.ar(
			// 	saturated, //
			// 	ironed, // ctr signal -> input, but pre-saturation
			// 	thresh: sat_threshold.clip(0.1, 1),
			// 	slopeBelow: 1,  // 1 means no comp pre-knee
			// 	slopeAbove: 0.5, // post-knee
			// 	clampTime: 0.01, // fast attack
			// 	relaxTime: 0.1 // fast release
			// );

			compressed = saturated;

		([
			phased, mixed/2,
			// phaseSlicedRm, phaseSlicedAm,
			phase/4,
			phaseRm,
			// crossingSliced, counterSliced/mod,
			// crossing, counter/mod,
			signal1, signal2, signal3,
			signal4
		].scope(name: "QuiltScope", bufsize: 4096*2));

			Out.ar(0, Pan2.ar(compressed, pan * (1 - (pan_lfo_amount * SinOsc.kr(pan_lfo_freq, pan_lfo_phase, 0.5, 0.5)))));
		}).add;


def.send(s);
s.sync;

~synth = Synth.new(\Quilt);
g_freq = ~updateScopeRange.(g_freq, g_mod1);

win = Window("Synth Controls", Rect(100, 100, 350, 450)).front;
win.onClose = {
	Stethoscope.ugenScopes.do({ arg scope, i; scope.quit() });
    ~synth.free;
	~sawBuffer.free;
	~sawValues.free;
};


// ------------------------------------
// controls - main

ly = 10;

StaticText(win, Rect(10, ly, 50, 20)).string_("f");
freqLabel = StaticText(win, Rect(270, ly, 50, 20)).string_(d_freq);
freqSlider = Slider(win, Rect(70, ly, 200, 20));
freqSlider.action = { |slider|
	g_freq = slider.value.linexp(0, 1, 20, 2000);
	g_freq = ~updateScopeRange.(g_freq, g_mod1);
	~synth.set(\freq, g_freq);
	freqLabel.string = g_freq;
};
freqSlider.valueAction = d_freq.explin(20, 2000, 0, 1);
ly = ly + lh;

ampLabel = StaticText(win, Rect(10, ly, 50, 20));
ampLabel.string = "a";
ampSlider = Slider(win, Rect(70, ly, 200, 20));
ampSlider.action = { |slider|
    ~synth.set(\amp_offset, slider.value);
};
ampSlider.valueAction = d_amp_offset;
ly = ly + lh;


StaticText(win, Rect(10, ly, 50, 20)).string_("cutoff");
Slider(win, Rect(70, ly, 200, 20))
.action_ ( { |slider|
	~synth.set(\phased_cutoff, slider.value.linexp(0, 1, 20, 20000));
	})
.valueAction_(15000.explin(20, 20000, 0, 1));
ly = ly + lh;

// ------------------------------------
// controls - modulators

StaticText(win, Rect(10, ly, 50, 20)).string_("pm f");
Slider(win, Rect(70, ly, 200, 20))
.action_ ( { |slider|
	~synth.set(\pmFreq, slider.value.linexp(0, 1, 0.1, 20000));
	})
.valueAction_(0);
ly = ly + lh;

StaticText(win, Rect(10, ly, 50, 20)).string_("pm a");
Slider(win, Rect(70, ly, 200, 20))
.action_ ( { |slider|
	~synth.set(\pmAmt, slider.value);
	})
.valueAction_(0);
ly = ly + lh;

StaticText(win, Rect(10, ly, 50, 20)).string_("m1 p");
Slider(win, Rect(70, ly, 200, 20))
.action_ ( { |slider|
	var p = slider.value.linlin(0, 1, -1, 1);
	~synth.set(\syncPhase, p);
	})
.valueAction_(0.5);
ly = ly + lh;

StaticText(win, Rect(10, ly, 50, 20)).string_("m2 p");
Slider(win, Rect(70, ly, 200, 20))
.action_ ( { |slider|
	var p = slider.value.linlin(0, 1, -1, 1);
	~synth.set(\syncPhaseSliced, p);
	})
.valueAction_(0.5);
ly = ly + lh;

//mod1Label = StaticText(win, Rect(10, ly, 50, 20));
//mod1Label.string = "m1";
StaticText(win, Rect(10, ly, 50, 20)).string_("m1");
mod1Slider = Slider(win, Rect(70, ly, 200, 20));
mod1Slider.action = { |slider|
	g_mod1 = slider.value.linlin(0, 1, 2, 15).round;
	g_freq = ~updateScopeRange.(g_freq, g_mod1);
	~synth.set(\freq, g_freq);
    ~synth.set(\mod, g_mod1);
};
mod1Slider.valueAction = d_mod1.linlin(2, 15, 0, 1);
ly = ly + lh;

mod1AmtLabel = StaticText(win, Rect(10, ly, 50, 20));
mod1AmtLabel.string = "m1 a";
mod1AmtSlider = Slider(win, Rect(70, ly, 200, 20));
mod1AmtSlider.action = { |slider|
    ~synth.set(\npolarProj, slider.value);
};
mod1AmtSlider.valueAction = d_mod1a;
ly = ly + lh;

mod1FreqLabel = StaticText(win, Rect(10, ly, 50, 20));
mod1FreqLabel.string = "m1 f";
mod1FreqSlider = Slider(win, Rect(70, ly, 200, 20));
mod1FreqSlider.action = { |slider|
    ~synth.set(\npolarRotFreq, slider.value.linexp(0, 1, 0.1, 20000));
};
mod1FreqSlider.valueAction = d_mod1f.explin(0.1, 20000, 0, 1);
ly = ly + lh;

mod2Label = StaticText(win, Rect(10, ly, 50, 20));
mod2Label.string = "m2";
mod2Slider = Slider(win, Rect(70, ly, 200, 20));
mod2Slider.action = { |slider|
    ~synth.set(\syncRatio, slider.value.linlin(0, 1, 1, 8).round);
};
mod2Slider.valueAction = d_mod2.linlin(1, 8, 0, 1);
ly = ly + lh;

mod2AmtLabel = StaticText(win, Rect(10, ly, 50, 20));
mod2AmtLabel.string = "m2 a";
mod2AmtSlider = Slider(win, Rect(70, ly, 200, 20));
mod2AmtSlider.action = { |slider|
    ~synth.set(\npolarProjSliced, slider.value);
};
mod2AmtSlider.valueAction = d_mod2a;
ly = ly + lh;

mod2FreqLabel = StaticText(win, Rect(10, ly, 50, 20));
mod2FreqLabel.string = "m2 f";
mod2FreqSlider = Slider(win, Rect(70, ly, 200, 20));
mod2FreqSlider.action = { |slider|
    ~synth.set(\npolarRotFreqSliced, slider.value.linexp(0, 1, 0.1, 20000));
};
mod2FreqSlider.valueAction = d_mod2f.explin(0.1, 20000, 0, 1);
ly = ly + lh;


// ------------------------------------
// controls - indices

i1Label = StaticText(win, Rect(10, ly, 50, 20));
i1Label.string = "i1";
i1Slider = Slider(win, Rect(70, ly, 200, 20));
i1Slider.action = { |slider|
    ~synth.set(\index1, slider.value.linlin(0, 1, 0, 3).round);
};
ly = ly + lh;

i2Label = StaticText(win, Rect(10, ly, 50, 20));
i2Label.string = "i2";
i2Slider = Slider(win, Rect(70, ly, 200, 20));
i2Slider.action = { |slider|
    ~synth.set(\index2, slider.value.linlin(0, 1, 0, 3).round);
};
ly = ly + lh;

i3Label = StaticText(win, Rect(10, ly, 50, 20));
i3Label.string = "i3";
i3Slider = Slider(win, Rect(70, ly, 200, 20));
i3Slider.action = { |slider|
    ~synth.set(\index3, slider.value.linlin(0, 1, 0, 3).round);
};
ly = ly + lh;

i4Label = StaticText(win, Rect(10, ly, 50, 20));
i4Label.string = "i4";
i4Slider = Slider(win, Rect(70, ly, 200, 20));
i4Slider.action = { |slider|
    ~synth.set(\index4, slider.value.linlin(0, 1, 0, 3).round);
};
ly = ly + lh;

};
)
