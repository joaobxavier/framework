package nl.tudelft.bt.model.work.carlos;

import java.awt.Color;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class LiveDeadBiomassSpecies extends BiomassSpecies {
	protected SoluteSpecies _autoInducer;

	protected double _threshold;


	public LiveDeadBiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume, SoluteSpecies autoInducer,
			double threshold) throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
		_threshold = threshold;
		_autoInducer = autoInducer;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.BiomassSpecies#createBiomassParticle()
	 */
	public BiomassParticle createBiomassParticle() {
		return new LiveDeadBiomassParticle(this);
	}
}
