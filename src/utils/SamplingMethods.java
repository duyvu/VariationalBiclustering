package utils;

import umontreal.iro.lecuyer.rng.RandomStream;

public class SamplingMethods {

	/* Draw a sample from a Bernoulli distribution p */
	public static int drawBernoulli(RandomStream sampler, double p) {
		return (sampler.nextDouble() > p) ? 0 : 1;
	}

	/* Draw a single sample from a multinomial distribution p */
	public static int drawSingleMultinomial(RandomStream sampler, double[] p) {
		double b = 0, r = sampler.nextDouble();
		for (int i = 0; i < p.length; i++) {
			b += p[i];
			if (b > r) {
				return i;
			}
		}
		return p.length - 1;
	}

	/* Draw a sample from a multinomial distribution n, p */
	public static int[] drawMultipleMultinomial(RandomStream sampler, int n,
			double[] p) {
		int[] sample = new int[p.length];
		for (int k = 0; k < p.length; k++)
			sample[k] = 0;
		for (int i = 0; i < n; i++)
			sample[drawSingleMultinomial(sampler, p)]++;
		return sample;
	}

}
