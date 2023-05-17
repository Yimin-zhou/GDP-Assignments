package workshop;

import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.List;
import java.awt.Panel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Vector;

import Jama.Matrix;
import Jama.SingularValueDecomposition;
import jv.geom.PgElementSet;
import jv.number.PuDouble;
import jv.number.PuInteger;
import jv.object.PsConfig;
import jv.object.PsDialog;
import jv.object.PsUpdateIf;
import jv.objectGui.PsList;
import jv.project.PgGeometryIf;
import jv.project.PvGeometryIf;
import jv.viewer.PvDisplay;
import jvx.project.PjWorkshop_IP;

import jv.vecmath.PdVector;
import jv.vecmath.PdMatrix;
import java.util.Random;



/**
 * Info Panel of Workshop for surface registration
 *
 */
public class Registration_IP extends PjWorkshop_IP implements ActionListener{
	protected	List			m_listActive;
	protected	List			m_listPassive;
	protected	Vector			m_geomList;
	protected	Registration	m_registration;
	protected   Button			m_bSetSurfaces;

	protected 	PuDouble 		m_threshold;

	protected 	PuInteger		m_icp_points;

	protected 	PuInteger		m_icp_iteration;

	protected Button m_bStartICP;
	// task 2.1 step 1
	private int icp_points = 50;
	private Random random = new Random();
	private 	double 		k = 1.5;
	private int icp_iteration = 10;


	/** Constructor */
	public Registration_IP () {
		super();
		if (getClass() == Registration_IP.class)
			init();
	}

	/**
	 * Informational text on the usage of the dialog.
	 * This notice will be displayed if this info panel is shown in a dialog.
	 * The text is split at line breaks into individual lines on the dialog.
	 */
	public String getNotice() {
		return "This text should explain what the workshop is about and how to use it.";
	}

	/** Assign a parent object. */
	public void setParent(PsUpdateIf parent) {
		super.setParent(parent);
		m_registration = (Registration)parent;

		addSubTitle("Select Surfaces to be Registered");

		Panel pGeometries = new Panel();
		pGeometries.setLayout(new GridLayout(1, 2));

		Panel Passive = new Panel();
		Passive.setLayout(new BorderLayout());
		Panel Active = new Panel();
		Active.setLayout(new BorderLayout());
		Label ActiveLabel = new Label("Surface P");
		Active.add(ActiveLabel, BorderLayout.NORTH);
		m_listActive = new PsList(5, true);
		Active.add(m_listActive, BorderLayout.CENTER);
		pGeometries.add(Active);
		Label PassiveLabel = new Label("Surface Q");
		Passive.add(PassiveLabel, BorderLayout.NORTH);
		m_listPassive = new PsList(5, true);
		Passive.add(m_listPassive, BorderLayout.CENTER);
		pGeometries.add(Passive);
		add(pGeometries);

		Panel pSetSurfaces = new Panel(new BorderLayout());
		m_bSetSurfaces = new Button("Set selected surfaces");
		m_bSetSurfaces.addActionListener(this);
		pSetSurfaces.add(m_bSetSurfaces, BorderLayout.CENTER);
		add(pSetSurfaces);

		// task 2
		m_threshold = new PuDouble("Threshold");
		m_threshold.setDefBounds(-50,50,0.1,1);
		m_threshold.addUpdateListener(this);
		m_threshold.init();
		add(m_threshold.getInfoPanel());

		m_icp_points = new PuInteger("ICP Points");
		m_icp_points.setDefBounds(1,100,1,1);
		m_icp_points.addUpdateListener(this);
		m_icp_points.init();
		add(m_icp_points.getInfoPanel());

		m_icp_iteration = new PuInteger("ICP Iteration");
		m_icp_iteration.setDefBounds(1,50,1,1);
		m_icp_iteration.addUpdateListener(this);
		m_icp_iteration.init();
		add(m_icp_iteration.getInfoPanel());

		m_bStartICP = new Button("Start ICP");
		m_bStartICP.addActionListener(this);
		Panel pStartICP = new Panel(new BorderLayout());
		pStartICP.add(m_bStartICP, BorderLayout.CENTER);
		add(pStartICP);

		updateGeomList();
		validate();
	}

	/** Initialisation */
	public void init() {
		super.init();
		setTitle("Surface Registration");

	}

	/** Set the list of geometries in the lists to the current state of the display. */
	public void updateGeomList() {
		Vector displays = m_registration.getGeometry().getDisplayList();
		int numDisplays = displays.size();
		m_geomList = new Vector();
		for (int i=0; i<numDisplays; i++) {
			PvDisplay disp =((PvDisplay)displays.elementAt(i));
			PgGeometryIf[] geomList = disp.getGeometries();
			int numGeom = geomList.length;
			for (int j=0; j<numGeom; j++) {
				if (!m_geomList.contains(geomList[j])) {
					//Take just PgElementSets from the list.
					if (geomList[j].getType() == PvGeometryIf.GEOM_ELEMENT_SET)
						m_geomList.addElement(geomList[j]);
				}
			}
		}
		int nog = m_geomList.size();
		m_listActive.removeAll();
		m_listPassive.removeAll();
		for (int i=0; i<nog; i++) {
			String name = ((PgGeometryIf)m_geomList.elementAt(i)).getName();
			m_listPassive.add(name);
			m_listActive.add(name);
		}
	}

	public boolean update(Object event) {
		if (event == m_threshold) {
			k = m_threshold.getValue();
			return true;
		} else if (event == m_icp_points) {
			icp_points = m_icp_points.getValue();
			return true;
		} else if (event == m_icp_iteration) {
			icp_iteration = m_icp_iteration.getValue();
			return true;
		}else
			return super.update(event);
	}
	/**
	 * Handle action events fired by buttons etc.
	 */
	public void actionPerformed(ActionEvent event) {
		Object source = event.getSource();
		if (source == m_bSetSurfaces) {
			m_registration.setGeometries((PgElementSet)m_geomList.elementAt(m_listActive.getSelectedIndex()),
					(PgElementSet)m_geomList.elementAt(m_listPassive.getSelectedIndex()));
			return;
		}
		else if (source == m_bStartICP) {
			performICP();
			return;
		}
	}

	// task 2
	public void performICP() {
		PgElementSet p = (PgElementSet) m_geomList.elementAt(m_listActive.getSelectedIndex());
		PgElementSet q = (PgElementSet) m_geomList.elementAt(m_listPassive.getSelectedIndex());

		for (int z = 0; z < icp_iteration; z++)
		{
			int[] selectedIndices = selectRandomVertices(p, icp_points);
			ArrayList<PdVector> pVertices = new ArrayList<>();
			ArrayList<PdVector> qVertices = new ArrayList<>();
			ArrayList<Double> distances = new ArrayList<>();

			for (int i : selectedIndices) {
				PdVector pi = p.getVertex(i);
				PdVector qi = findClosestVertex(pi, q);

				// Store vertices and distance
				pVertices.add(pi);
				qVertices.add(qi);
				distances.add(pi.dist(qi));
			}

			Collections.sort(distances);
			double medianDistance;
			if (distances.size() % 2 == 0) {
				medianDistance = (distances.get(distances.size() / 2 - 1) + distances.get(distances.size() / 2)) / 2.0;
			} else {
				medianDistance = distances.get(distances.size() / 2);
			}

			double threshold = k * medianDistance;

			System.out.println("k: " + k);
			System.out.println("median distance: " + medianDistance);

			ArrayList<PdVector> newPVertices = new ArrayList<>();
			ArrayList<PdVector> newQVertices = new ArrayList<>();

			for (int i = 0; i < distances.size(); i++) {
				if (distances.get(i) <= threshold) {
					newPVertices.add(pVertices.get(i));
					newQVertices.add(qVertices.get(i));
				}
			}

			pVertices = newPVertices;
			qVertices = newQVertices;

			// task 2.3
			// Compute centroids of P and Q points
			PdVector centroidP = new PdVector(3);
			PdVector centroidQ = new PdVector(3);
			for (int i = 0; i < pVertices.size(); i++) {
				centroidP.add(pVertices.get(i));
				centroidQ.add(qVertices.get(i));
			}
			centroidP.multScalar(1.0 / pVertices.size());
			centroidQ.multScalar(1.0 / qVertices.size());

			// Build matrix M
			PdMatrix M = new PdMatrix(3, 3);
			for (int i = 0; i < pVertices.size(); i++) {
				PdVector pCentroid = PdVector.subNew(pVertices.get(i), centroidP);
				PdVector qCentroid = PdVector.subNew(qVertices.get(i), centroidQ);
				PdMatrix pqT = outerProductNew(pCentroid, qCentroid);
				M.add(pqT);
			}

			// compute SVD of M
			Matrix jamaM = new Matrix(M.getEntries());
			SingularValueDecomposition SVD = jamaM.svd();

			// compute optimal rotation matrix R
			// compute the determinant of VU^T
			Matrix VUt = SVD.getV().times(SVD.getU().transpose());
			double det = VUt.det();

			Matrix correction = Matrix.identity(3, 3);
			correction.set(2, 2, det);

			// compute optimal rotation matrix R
			Matrix R_jama = SVD.getV().times(correction).times(SVD.getU().transpose());
			// convert Jama matrix to PdMatrix
			PdMatrix R = new PdMatrix(R_jama.getArray());

			// compute translation vector t
			PdVector t = new PdVector(3);
			t = PdVector.subNew(centroidQ, R.leftMultMatrix(t, centroidP));

			// Apply the optimal rigid transformation to P
			for (int i = 0; i < p.getNumVertices(); i++) {
				PdVector pi = p.getVertex(i);
				PdVector qi = R.leftMultMatrix(pi, pi);
				qi.add(t);
				p.setVertex(i, qi);
			}

			// update the display
			p.update(p);
		}
	}

	private int[] selectRandomVertices(PgElementSet mesh, int numVertices) {
		int[] indices = new int[numVertices];
		for (int i = 0; i < numVertices; i++) {
			indices[i] = random.nextInt(mesh.getNumVertices());
		}
		return indices;
	}

	private PdVector findClosestVertex(PdVector point, PgElementSet mesh) {
		PdVector closest = null;
		double closestDistance = Double.MAX_VALUE;

		// brute force search
		for (int i = 0; i < mesh.getNumVertices(); i++) {
			PdVector vertex = mesh.getVertex(i);
			double distance = point.dist(vertex);
			if (distance < closestDistance) {
				closest = vertex;
				closestDistance = distance;
			}
		}

		return closest;
	}
	public static PdMatrix outerProductNew(PdVector a, PdVector b) {
		PdMatrix result = new PdMatrix(3, 3);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				result.setEntry(i, j, a.getEntry(i) * b.getEntry(j));
			}
		}
		return result;
	}

	/**
	 * Get information which bottom buttons a dialog should create
	 * when showing this info panel.
	 */
	protected int getDialogButtons()		{
		return PsDialog.BUTTON_OK;
	}
}