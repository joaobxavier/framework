package nl.tudelft.bt.model.work.carlos;

import java.awt.Color;
import java.util.Iterator;
import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Height2MassDetachment;
import nl.tudelft.bt.model.detachment.levelset.functions.Radius2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BiofilmBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.FixedOnTopBoundaryLayer;
import nl.tudelft.bt.model.multigrid.boundary_layers.SphericalDilationBoundaryLayer;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.particlebased.tumor.TumorModelHandler;
import nl.tudelft.bt.model.reaction.*;

/**
 * @author Carlos Carmona Fontaine (carmonac@mskcc.org)
 */
public class SpatialStructureGradientsNGSandwich extends ModelHandler {
	// public class Example1KWithTumorHandler extends TumorModelHandler {
	// All the model parameters are defined here as static attributes
	// at the beginning of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	protected static String outputDirectory = "/Users/xavierj/results/sandwich";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.
	// The output directory is were the program will store all the results.
	// Choose a path to an existing folder in your system.
	// EXAMPLE: if you choose "e:/results/example1/" directory "e:\results" must
	// exist in your computer. The subdirectory "example1" will be created
	// if it is non-existant. If it exists, its contents will be deleted
	// during the program initialization

	// geometry (default is 2D) - change value to 3 for 3D
	protected static int geometry = 2;

	/**
	 * Constants
	 */

	// /////////// Computation parameters

	protected static float systemSize = 500; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 5f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.6f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 0f;

	// other model parameters
	protected static int gridSide = 33; // multigrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less] NEED
											// THIS?

	protected static float rdetach = 1e-10f; // NEED THIS?

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 1000;

	// /////////// Biochemical parameters //

	protected static SoluteSpecies oxygen;
	protected static SoluteSpecies lactate;

	protected static float oxygenBulkConcentration = 8e-3f; // [g/L]Joao's paper 
															// Should be 16-22 ml O2/100ml arterial blood
															// and 13-16 ml in venous blood ()
	protected static float lactateBulkConcentration = 0; // [g/L]

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	// Glucose diffusivity (Casciari et al 1988) 1.1e-6 cm2/s = 3.96e5um2/h
	// mouse; 5.5e-7-2.3e-7 cm2/s = 1.98e4um2/h to 8.2e3 um2/h
	
	//Use this line to turn off gradients
	//private static float oxygenDiffusivity = 7.2e7f; // 10x Diffusivity

	//Use this line to turn on gradients
	private static float oxygenDiffusivity = 7.2e6f;//[um^2/h]Thomlinson: 2e-5
														// cm2/s
	private static float lactateDiffusivity = 1.04e6f; // [um^2/h] 2.9e-6 cm2/s.
														// Measurement of
														// Effective
														// Diffusivities of
														// Lactose and Lactic
														// Acid in 3% Agarose
														// Gel Membrane
														// Amarjeet S. Bassi,
														// Sohrab Rohani and D.
														// G. Macdonaid

	// protected static SoluteSpecies Hydrogen;

	/**
	 * Variables
	 */

	// /////////// Biochemical parameters

	// Yield coefficients
	//private static float gYo2 = 0.8f; // [gX/gS] Glucose Yield Aerobic
	private static float gYana = 0.5f; // [gX/gS] Lactate Yield 
	private static float oY = 0.5f; // [gX/gS] Oxygen Yield
	private static float lY = 1f - gYana; // [gX/gS] Lactate production Yield

	// Satutration constants

	private static float oKST = oxygenBulkConcentration / 2f; // [g/L]
	private static float oKSM = oxygenBulkConcentration / 20f; // [g/L]

	private static float lKST = 1e-1f; // [g/L] saturation constant lactate
										// effect on tumor cells
	private static float lKSM = 1e-1f; // [g/L] saturation constant lactate
										// effect on macrophages

	// /////////// Cell parameters

	protected static int numberofCellGroups = 2;
	//
	// Particulate species (biomass X)
	protected static float specificMassT = 100f; // [g/L]
	protected static float specificMassM = 100f; // [g/L]

	// Max growth rates
	protected static float uMaxT = 0.01f; // [1/h]
	protected static float uMaxM = uMaxT * 1f; // [1/h]

	// Death rates
	protected static float kdT = uMaxT * 4f;// lactate effect strength on tumor cells (relative to growth)
	protected static float kdM = kdT * 2f;// lactate effect strength on macrophages (relative to cancer sensitivity)

	
	/**
	 * Definitions and Processes
	 */

	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// 1. Create the solutes
		// Oxygen
		oxygen = new SoluteSpecies("Oxygen", oxygenDiffusivity);
		// set up the simplest type of bulk concentration: constant
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		// Glucose

		// Lactate
		lactate = new SoluteSpecies("Lactate", lactateDiffusivity);
		// set up the simplest type of bulk concentration: constant
		lactate.setBulkConcentration(new ConstantBulkConcentration(
				lactateBulkConcentration));

		// 2. Create the particulate species (solids)

		// X active mass
		ParticulateSpecies activeT = new ParticulateSpecies("activeT",
				specificMassT, Color.blue);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spT = { activeT };
		float[] fractionalVolumeCompositionT = { 1.0f };

		ParticulateSpecies activeM = new ParticulateSpecies("activeM",
				specificMassM, Color.green);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spM = { activeM };
		float[] fractionalVolumeCompositionM = { 1.0f };

		// 3. Create the biomass species
		BiomassSpecies speciesT = new BiomassSpecies("speciesT", spT,
				fractionalVolumeCompositionT);
		speciesT.setActiveMass(activeT);
		// speciesCC.getColorFromGrowth();

		BiomassSpecies speciesM = new BiomassSpecies("speciesM", spM,
				fractionalVolumeCompositionM);
		speciesM.setActiveMass(activeM);
		// speciesTAM.getColorFromGrowth();

		// 4. Create the Reaction factors, Monod and inhibition coefficients

		// Growth factor production

		// Reaction egfsecretion = new Reaction("egfsecretion", activeTAM, iegf,
		// 0); //species or active?

		// ProcessFactor mSo = new Saturation(glucose, oKS);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate
		// growth
		// aerobic growth weight
		ProcessFactor TumorAeroWeight = new Saturation(oxygen, oKST);
		// anaerobic growth weight
		ProcessFactor TumorAnaWeight = new Inhibition(oxygen, oKST);

		// aerobic growth weight
		ProcessFactor MacroAeroWeight = new Saturation(oxygen, oKSM);
		// anaerobic growth weight
		ProcessFactor MacroAnaWeight = new Inhibition(oxygen, oKSM);

		// ProcessFactor Endogenous = new Saturation(lactate, lKST);
		ProcessFactor LactateEffT = new Saturation(lactate, lKST);
		ProcessFactor LactateEffM = new Saturation(lactate, lKSM);

		// 5. Create the reactions
		// Aerobic growth
		Reaction aeroGrowthT = new Reaction("aeroGrowthT", activeT, uMaxT, 1);
		aeroGrowthT.addFactor(TumorAeroWeight);

		Reaction aeroGrowthM = new Reaction("aeroGrowthM", activeM, uMaxM, 1);
		aeroGrowthM.addFactor(MacroAeroWeight);

		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass

		// Anaerobic growth
		Reaction anaGrowthT = new Reaction("anaGrowthT", activeT, uMaxT, 1);
		anaGrowthT.addFactor(TumorAnaWeight);

		Reaction anaGrowthM = new Reaction("anaGrowthM", activeM, uMaxM, 1);
		anaGrowthM.addFactor(MacroAnaWeight);

		// Endogenous Decay and death

		Reaction LactateKillingT = new Reaction("LactateKillingT", activeT,
				kdT, 1);
		LactateKillingT.addFactor(LactateEffT);

		Reaction LactateKillingM = new Reaction("LactateKillingM", activeM,
				kdM, 1);
		LactateKillingM.addFactor(LactateEffM);

		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactiveT = new NetReaction(3);
		rsXactiveT.addReaction(aeroGrowthT, 1);
		rsXactiveT.addReaction(anaGrowthT, 1);
		rsXactiveT.addReaction(LactateKillingT, -1);
		activeT.setProcesses(rsXactiveT);

		NetReaction rsXactiveM = new NetReaction(3);
		rsXactiveM.addReaction(aeroGrowthM, 1);
		rsXactiveM.addReaction(anaGrowthM, 1);
		rsXactiveM.addReaction(LactateKillingM, -1);
		activeM.setProcesses(rsXactiveM);
		// This defines that biomass growth rate is 1*rX
		//
		// assign reaction stoichiometry to the solutes
		// substrate

		NetReaction rsOxygen = new NetReaction(2);
		rsOxygen.addReaction(aeroGrowthT, -(1 / oY));
		rsOxygen.addReaction(aeroGrowthM, -(1 / oY));
		oxygen.setProcesses(rsOxygen);

		NetReaction rsLactate = new NetReaction(2);
		rsLactate.addReaction(anaGrowthT, (lY));
		rsLactate.addReaction(anaGrowthM, (lY));
		lactate.setProcesses(rsLactate);

		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesT);
		addBiomassSpecies(speciesM);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(lactate);

	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
	}

	/*
	 * (non-Javadoc)
	 */
	protected void inoculate() {
		int[] nCells = { initialParticleNumber, initialParticleNumber };
		// inoculateRandomlyInsideRadius(nCells, initialColonyRadius);
	//	inoculateRandomlyEverywhere(nCells);
		inoculateRandomly(nCells);

	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Height2MassDetachment(0);
		setDetachmentHandler(df);
	}

	protected void createBoundaryLayer(float h)
			throws MultigridSystemNotSetException {
		_boundaryLayer = new FixedOnTopBoundaryLayer();
		// create the boundary conditions
		MultigridVariable
				.setBoundaryConditions(new BiofilmBoundaryConditions());
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new SpatialStructureGradientsNGSandwich();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// The following code will be omitted if no vizuals are desired
		// start decorationg the application
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		app = new SeriesVizualizer(app, biomass);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// --- nothing to set in this case: constant bulk concentration
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			app.initializeDiffusionReactionSystem(); // also innoculates

			// grafics
			app.addStateWriter(new PovRayWriter());
			app.addStateWriter(new ImageWriter());
			app.addStateWriter(new ImageWriter(lactate));
			app.addStateWriter(new ImageWriter(oxygen));
			// app.addStateWriter(new ImageWriter(egf));
			// app.addStateWriter(new ImageWriter(csf1));
			// app.addStateWriter(new ImageWriter(IL4));
			app.addStateWriter(new SoluteConcentrationWriter());
			app.addStateWriter(new SolidsConcentrationWriter());
			app.addStateWriter(new ParticlePositionWriter());
			// app.addStateWriter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getErodedBiomassSeries());
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getSloughedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
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
			System.exit(-1);
		}
		try {
			// start iterating cycle
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		System.out.println("Simulation finished.");
	}
}