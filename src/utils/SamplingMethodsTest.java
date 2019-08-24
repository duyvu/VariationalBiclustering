package utils;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import umontreal.iro.lecuyer.rng.LFSR113;
import umontreal.iro.lecuyer.rng.RandomStream;

public class SamplingMethodsTest {

	protected RandomStream sampler;
	protected double[] p;

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		sampler = new LFSR113();
		p = new double[4];
		p[0] = .1;
		p[1] = .2;
		p[2] = .3;
		p[3] = .4;
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testDrawSingleMultinomial() {
		for (int i = 0; i < 30; i++)
			System.out.println(SamplingMethods
					.drawSingleMultinomial(sampler, p));
	}

	@Test
	public void testDrawMultipleMultinomial() {
		for (int i = 0; i < 30; i++) {
			int[] sample = SamplingMethods.drawMultipleMultinomial(sampler, 10,
					p);
			for (int k = 0; k < p.length; k++)
				System.out.print(sample[k] + "\t");
			System.out.println();
		}
	}

}
