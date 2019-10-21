package nl.tudelft.bt.model.work.drift;

import java.awt.Color;
import java.io.IOException;
import java.util.Iterator;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Radius2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.reaction.*;

/**
 * Simulates growth of two strains with neutral advantage in a 2D colony
 * spreading on agar plate. Nutrients diffuse into colony from outside with no
 * external mass transfer resistance. Physiology is implemented with a threshold
 * concentration at S = 0.
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu)
 */
public class DriftIn2DColony extends GranuleModelHandler {
	// All the model parameters are defined here as static attributes
	// at the beginning of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory = "C:/results/colony/drift/test/";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.

	// Growth
	protected static float uMax = 1f; // [h-1]
	private static float YS = 1.5f; // [gCOD-PHB/gCOD-S]

	// geometry Must be 2 for 2D colony growth
	protected static int geometry = 2;

	// Nutrient properties
	protected static float sBulk = 3e-3f; // [gO2/L] DO 20
	private static float diffusivity = 4e5f; // [um2/h]

	// Particulate species (biomass) properties
	protected static float density = 150f; // [gCOD-H/L]

	// Computation parameters
	// Size of computational volume (size of size of square)
	protected static float systemSize = 500; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 1f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 10 um
	protected static float relativeBoundaryLayer = 0 / systemSize;

	protected static float maximumColonyRadius = 200; // [um]

	protected static float initialColonyRadius = 25; // [um]

	// other model parameters
	protected static int gridSide = 129; // multigrid grid side

	// Don't change this
	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// Don't change this, detachment rate, leave at zero to form round granules
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 100;

	// iteration finish time
	protected static float simulationFinishTime = 48f; // [h]

	// outpute (write results to file) every:
	protected static float outputEvery = 0.2f; // [h]

	// /END OF PARAMETERS

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	protected void defineSpeciesAndReactions() throws ModelException {
		// substrate (S)
		SoluteSpecies substrate = new SoluteSpecies("substrate", diffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(sBulk));
		// 2. Create the particulate species (solids)
		// The active mass of green strain
		ParticulateSpecies activeGreen = new ParticulateSpecies("activeH",
				density, Color.green);
		// The active mass of red strain
		ParticulateSpecies activeRed = new ParticulateSpecies("activeRed",
				density, Color.red);
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		// green strain
		ParticulateSpecies[] spGreen = { activeGreen };
		float[] fractionalVolumeCompositionGreen = { 1.0f };
		// 3. Create the biomass species
		BiomassSpecies speciesGreen = new BiomassSpecies("speciesGreen",
				spGreen, fractionalVolumeCompositionGreen);
		speciesGreen.setActiveMass(activeGreen);
		// red strain
		ParticulateSpecies[] spRed = { activeRed };
		float[] fractionalVolumeCompositionRed = { 1.0f };
		// 3. Create the biomass species
		BiomassSpecies speciesRed = new BiomassSpecies("speciesRed", spRed,
				fractionalVolumeCompositionRed);
		speciesRed.setActiveMass(activeRed);
		// 4. Create the Reaction factors
		ProcessFactor step = new Step(substrate, 0);
		// 5. Create the reactions
		// aerobic growth Green
		Reaction aerobicGrowthGreen = new Reaction("growthGreen", activeGreen,
				uMax, 1);
		aerobicGrowthGreen.addFactor(step);
		// aerobic growth Red
		Reaction aerobicGrowthRed = new Reaction("growthRed", activeRed, uMax,
				1);
		aerobicGrowthRed.addFactor(step);
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active Green
		NetReaction rsActiveGreen = new NetReaction(1);
		rsActiveGreen.addReaction(aerobicGrowthGreen, 1);
		activeGreen.setProcesses(rsActiveGreen);
		// active Red
		NetReaction rsActiveRed = new NetReaction(1);
		rsActiveRed.addReaction(aerobicGrowthRed, 1);
		activeRed.setProcesses(rsActiveRed);
		// assign reaction stoichiometry to the solutes
		// substrate (S)
		NetReaction rsSubstrate = new NetReaction(2);
		rsSubstrate.addReaction(aerobicGrowthGreen, -1 / YS);
		rsSubstrate.addReaction(aerobicGrowthRed, -1 / YS);
		substrate.setProcesses(rsSubstrate);
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesGreen);
		addBiomassSpecies(speciesRed);
		addSoluteSpecies(substrate);
	}

	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
	}

	/*
	 * (non-Javadoc)
	 */
	protected void inoculate() {
		int[] nCells = { initialParticleNumber, initialParticleNumber };
		inoculateRandomlyInsideRadius(nCells, initialColonyRadius);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Radius2MassDetachment(kdetach);
		setDetachmentHandler(df);
		// set the maximum granule radius
		try {
			Model.model().setVerticalCutoffSize(maximumColonyRadius);
		} catch (InvalidValueException e) {
			System.out.println(e);
			System.exit(-1);
		}
	}

	/**
	 * @param app
	 */
	protected static void setSystemParametersAndInitializeSystemSpace(
			ApplicationComponent app) {
		// create the space
		app.setSystemSpaceParameters(geometry, systemSize,
				relativeMaximumRadius, relativeMinimumRadius,
				relativeBoundaryLayer, gridSide, kShov);
		// initialize
		app.initializeSystemSpace();
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line
		if (args.length < 4) {
			throw new RuntimeException("input arguments missing:\n"
					+ "1: output directory (CAUTION!!!"
					+ " directory will be erased\n"
					+ "2: seed for random number generator\n"
					+ "3: flag for running with graphics (1 on, 0 off)"
					+ "4: diffusion-reaction length scale (e.g. 5)");
		}
		// parse inputs
		outputDirectory = args[0];
		int seed = Integer.parseInt(args[1]);
		Model.model().setSeed(seed);
		boolean runWithGraphics = (Integer.parseInt(args[2]) == 1);
		uMax *= (Float.parseFloat(args[3]) / 4.898979484f);
		// set numerics for multigrid
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new DriftIn2DColony();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the biovolume series
		VariableSeries biovolume = new BiovolumeSeries();
		VariableSeries[] runLengthSeries = { new RunLengthXSeries(),
				new RunLengthYSeries(), new RunLengthZSeries() };
		// The following code will be omitted if no vizuals are desired
		if (runWithGraphics) {
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// finally, the controller must be the last decorator to add
			// app = new VizualModelControler(app);
		}
		try {
			// create the space
			setSystemParametersAndInitializeSystemSpace(app);
			// initialize
			app.intializeStateWriters(outputDirectory);
			// Pov witer is added twice
			app
					.addStateWriter(new TimedStateWriterDecorator(
							new PovRayWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SoluteConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SolidsConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new ParticlePositionWriter()));
			// app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(biovolume);
			spw.addSeries(runLengthSeries[0]);
			spw.addSeries(runLengthSeries[1]);
			spw.addSeries(runLengthSeries[2]);
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
			// add the time constraints writer
			app.addStateWriter(new TimeConstraintsWriter());
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
			// add bulk concentrations of all solutes as variable series
			for (Iterator i = Model.model().getSoluteSpecies().iterator(); i
					.hasNext();) {
				SoluteSpecies s = (SoluteSpecies) i.next();
				spw.addSeries(s.getBulkConcentrationSeries());
				spw.addSeries(s.getRateTimeSeries());
			}
			// add particulate global masses to write
			for (Iterator i = Model.model().getParticulateSpecies().iterator(); i
					.hasNext();) {
				ParticulateSpecies s = (ParticulateSpecies) i.next();
				spw.addSeries(s.getTotalMassSeries());
			}
		} catch (ModelException e) {
			System.out.println(e);
			e.printStackTrace();
			System.exit(-1);
		}
		try {
			// start iterating cycle
			Model.model().setCompulsoryTimeStep(outputEvery);
			Model.model().setFinishIterationTime(simulationFinishTime);
			Model.model().setMaxRunLength(50f);
			// start the iteration
			app.writeState(); // write iteration 0
			app.startIterating();
		} catch (Exception e1) {
			try {
				app.forceWriteState();
			} catch (IOException e2) {
				System.err.println("Error serializing state:");
				System.err.println("");
				e2.printStackTrace();
			}
			System.err.println("");
			System.err.println("Program failed due to :");
			System.err.println("");
			e1.printStackTrace();
			System.out.println(e1);
		}
		System.out.println("Simulation finished.");
		System.exit(0);
	}
}