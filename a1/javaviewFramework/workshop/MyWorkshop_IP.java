package workshop;

import java.awt.Button;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.*;

import jv.geom.PgBndPolygon;
import jv.geom.PgEdgeStar;
import jv.number.PuDouble;
import jv.object.PsConfig;
import jv.object.PsDialog;
import jv.object.PsUpdateIf;
import jv.vecmath.PiVector;
import jvx.project.PjWorkshop_IP;

import jv.vecmath.PdVector;
import jv.geom.PgElementSet;

public class MyWorkshop_IP extends PjWorkshop_IP implements ActionListener {

	protected Button m_bMakeRandomElementColors;
	protected Button m_bMakeRandomVertexColors;
	protected PuDouble m_xOff;
	protected TextField	m_textField;

	// task 1.1 calculate genus
	protected Button m_bCalculateGenus;
	protected TextField m_genusResult;

	// task 1.2 calculate closed volume
	protected Button m_bCalculateVolume;
    protected TextField m_volumeResult;

	// task 1.3 calculate number of connected components
	protected Button m_bCalculateComponents;
    protected TextField m_componentsResult;
	
	MyWorkshop m_ws;
	
	public MyWorkshop_IP() {
		super();
		if(getClass() == MyWorkshop_IP.class)
			init();
	}
	
	public void init() {
		super.init();
		setTitle("Mesh and Surface Analysis");
	}
	
	public String getNotice() {
		return "This workshop can be used to calculate perform mesh and surface analysis to retrieve properties such as genus, volume and amount of components.";
	}
	
	public void setParent(PsUpdateIf parent) {
		super.setParent(parent);
		m_ws = (MyWorkshop)parent;
	
		addSubTitle("Properties");

//		m_bMakeRandomElementColors = new Button("Random Element Colors");
//		m_bMakeRandomElementColors.addActionListener(this);
//		m_bMakeRandomVertexColors = new Button("Random Vertex Colors");
//		m_bMakeRandomVertexColors.addActionListener(this);
		m_textField = addTextField("Number of Triangles:", 20);
		m_textField.setEditable(false);
		m_textField.setText(String.valueOf(m_ws.m_geom.getNumElements()));
//		Panel panel1 = new Panel(new FlowLayout(FlowLayout.CENTER));
//		panel1.add(m_bMakeRandomElementColors);
//		panel1.add(m_bMakeRandomVertexColors);
//		add(panel1);

		m_xOff = new PuDouble("X Offset");
		m_xOff.setDefBounds(-10,10,0.1,1);
		m_xOff.addUpdateListener(this);
		m_xOff.init();
		add(m_xOff.getInfoPanel());

		// task 1.1 calculate genus
        m_bCalculateGenus = new Button("Calculate Genus");
        m_bCalculateGenus.addActionListener(this);
        m_genusResult = new TextField("Genus: ", 20);
        m_genusResult.setEditable(false);
        Panel panel2 = new Panel(new FlowLayout(FlowLayout.CENTER));
        panel2.add(m_bCalculateGenus);
        panel2.add(m_genusResult);
        add(panel2);

		// task 1.2 calculate closed volume
        m_bCalculateVolume = new Button("Calculate Volume");
        m_bCalculateVolume.addActionListener(this);
        m_volumeResult = new TextField("Volume: ", 20);
        m_volumeResult.setEditable(false);
        Panel panel3 = new Panel(new FlowLayout(FlowLayout.CENTER));
        panel3.add(m_bCalculateVolume);
        panel3.add(m_volumeResult);
        add(panel3);

		// task 1.3 calculate number of connected components
		m_bCalculateComponents = new Button("Calculate Components");
        m_bCalculateComponents.addActionListener(this);
        m_componentsResult = new TextField("Components: ", 20);
        m_componentsResult.setEditable(false);
        Panel panel4 = new Panel(new FlowLayout(FlowLayout.CENTER));
        panel4.add(m_bCalculateComponents);
        panel4.add(m_componentsResult);
        add(panel4);

		
		validate();
	}
	
	
	public boolean update(Object event) {
		if (event == m_xOff) {
			m_ws.setXOff(m_xOff.getValue());
			m_ws.m_geom.update(m_ws.m_geom);
			return true;
		} else
			return super.update(event);
	}
	
	/**
	 * Handle action events fired by buttons etc.
	 */
	public void actionPerformed(ActionEvent event) {
		Object source = event.getSource();
		if (source == m_bMakeRandomElementColors) {
			m_ws.makeRandomElementColors();
			m_ws.m_geom.update(m_ws.m_geom);
			return;
		}
		else if (source == m_bMakeRandomVertexColors) {
			m_ws.makeRandomVertexColors();
			m_ws.m_geom.update(m_ws.m_geom);
			return;
		}
		// task 1.1 calculate genus	
		else if (source == m_bCalculateGenus) {
            calculateGenus();
            return;
        }
		// task 1.2 calculate closed volume
		else if (source == m_bCalculateVolume) {
            calculateVolume();
            return;
        }
		// task 1.3 calculate number of connected components
        else if (source == m_bCalculateComponents) {
            calculateComponents();
            return;
        }
	}
	// task 1.1 calculate genus
	private void calculateGenus() {
		int V = m_ws.m_geom.getNumVertices();
		int E = m_ws.m_geom.getNumEdges();
		int F = m_ws.m_geom.getNumElements();
		int genus = 0;

		// assume no boundary loops
		genus = (int) (1.0f - (float)((V - E + F) / 2.0f));
		if (genus < 0) {
			genus = 0;
		}
		m_genusResult.setText("Genus: " + genus);

		// print v e f
		// System.out.println("V: " + V + " E: " + E + " F: " + F);
	}

	// task 1.2 calculate closed volume
    private void calculateVolume() {
        double volume = 0.0;

		// if the mesh is not closed, return
		if (m_ws.m_geom.getNumBoundaryEdges() > 0) {
			m_volumeResult.setText("Volume: " + "The mesh is not closed");
			return;
		}

        for (int i = 0; i < m_ws.m_geom.getNumElements(); i++) {
			PiVector element = m_ws.m_geom.getElement(i);
            PdVector v0 = m_ws.m_geom.getVertex(element.getEntry(0));
            PdVector v1 = m_ws.m_geom.getVertex(element.getEntry(1));
            PdVector v2 = m_ws.m_geom.getVertex(element.getEntry(2));

			PdVector crossProduct = new PdVector();
			crossProduct = crossProduct.crossNew(v1, v2);
			volume += v0.dot(crossProduct) / 6.0;
        }

        // Display the result
        m_volumeResult.setText("Volume: " + volume);
    }

	// task 1.3 calculate number of connected components
    private void calculateComponents() {
        int numComponents = getNumberOfConnectedComponents(m_ws.m_geom);
        m_componentsResult.setText("Components: " + numComponents);
    }

	private int getNumberOfConnectedComponents(PgElementSet geom) {
		int num_components = 0;
		LinkedList<Integer> queue = new LinkedList<Integer>();
		int element_count = m_ws.m_geom.getNumElements();
		boolean[] visited = new boolean[element_count];

		for (int i = 0; i < element_count; i++) {
			if (!visited[i]) {
				visited[i] = true;
				queue.add(i);

				while (queue.size() != 0) {
					int element = queue.poll();
					PiVector neighbours = m_ws.m_geom.getNeighbours()[element];

					for (int j = 0; j < neighbours.getSize(); j++) {
						int entry = neighbours.getEntry(j);
						if (entry < 0)
							continue;

						if (!visited[entry]) {
							visited[entry] = true;
							queue.add(entry);
						}
					}
				}
				num_components++;
			}
		}

		return num_components;
	}	

	/**
	 * Get information which bottom buttons a dialog should create
	 * when showing this info panel.
	 */
	protected int getDialogButtons()		{
		return PsDialog.BUTTON_OK;
	}
}
