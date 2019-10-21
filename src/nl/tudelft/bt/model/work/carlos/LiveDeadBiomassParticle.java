package nl.tudelft.bt.model.work.carlos;

import java.awt.Color;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

/**
 * Implements the mutator particle
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Oct 18, 2006
 */
public class LiveDeadBiomassParticle extends BiomassParticle {

	public LiveDeadBiomassParticle(LiveDeadBiomassSpecies s) {
		super(s);
	}

	// TODO check if these methods can be removed
	@Override
	public BiomassParticle divide() {
		LiveDeadBiomassParticle daughter = (LiveDeadBiomassParticle) super
				.divide();
		return daughter;
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}

	@Override
	public Color getColorCore() {
		// The problem was that setCurrentLocation was calle after getValue!
		MultigridVariable.setCurrentLocation(getCenter());
		double a = ((LiveDeadBiomassSpecies) _biomassSpecies)._autoInducer
				.getValue();
		double t = ((LiveDeadBiomassSpecies) _biomassSpecies)._threshold;

		float f = (float) (1 / (1 + Math.pow(t / a, 10)));


		return new Color(f, 1 - f, 0);
	}
}
