package workshop;

import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.List;
import java.awt.Panel;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.CheckboxGroup;
import java.awt.Checkbox;
import java.util.Vector;

import javax.swing.JFrame;

import jv.geom.PgElementSet;
import jv.object.PsConfig;
import jv.object.PsDialog;
import jv.object.PsUpdateIf;
import jv.objectGui.PsList;
import jv.project.PgGeometryIf;
import jv.project.PvGeometryIf;
import jv.viewer.PvDisplay;
import jvx.project.PjWorkshop_IP;
import jv.vecmath.PdMatrix;

import jv.number.PuDouble;


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

	protected   Button			m_itClosestPointButton;
	protected	CheckboxGroup 	m_pointOrPlaneRadioButton;
	protected	Checkbox 		m_point;
	protected	Checkbox 		m_plane;

	protected 	Panel 			m_matrixPanel;
	protected 	TextField		m_r1c1;
	protected 	TextField		m_r1c2;
	protected 	TextField		m_r1c3;
	protected 	TextField		m_r2c1;
	protected 	TextField		m_r2c2;
	protected 	TextField		m_r2c3;
	protected 	TextField		m_r3c1;
	protected 	TextField		m_r3c2;
	protected 	TextField		m_r3c3;
	protected   Button			m_apply;
	protected   Button			m_reset;
	protected   Button			m_applyRotation;
	protected 	PuDouble 		m_degrees;
	protected	CheckboxGroup 	m_axisRadioButton;
	protected	Checkbox 		m_yAxis;
	protected	Checkbox 		m_xAxis;
	protected	Checkbox 		m_zAxis;


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

		m_itClosestPointButton = new Button("Iterative closest point");
		m_itClosestPointButton.addActionListener(this);
		add(m_itClosestPointButton);
		m_pointOrPlaneRadioButton = new CheckboxGroup();
		m_point = new Checkbox("Point to point", m_pointOrPlaneRadioButton, true);
        m_plane = new Checkbox("Point to plane", m_pointOrPlaneRadioButton, false);
        add(m_point);
        add(m_plane);

        m_matrixPanel = new Panel();
      	m_matrixPanel.setLayout(new GridLayout(3,3));
      	m_r1c1 = new TextField("1");
		m_r1c2 = new TextField("0");
		m_r1c3 = new TextField("0");
		m_r2c1 = new TextField("0");
		m_r2c2 = new TextField("1");
		m_r2c3 = new TextField("0");
		m_r3c1 = new TextField("0");
		m_r3c2 = new TextField("0");
		m_r3c3 = new TextField("1");
		m_matrixPanel.add(m_r1c1);
		m_matrixPanel.add(m_r1c2);
		m_matrixPanel.add(m_r1c3);
		m_matrixPanel.add(m_r2c1);
		m_matrixPanel.add(m_r2c2);
		m_matrixPanel.add(m_r2c3);
		m_matrixPanel.add(m_r3c1);
		m_matrixPanel.add(m_r3c2);
		m_matrixPanel.add(m_r3c3);

      	add(m_matrixPanel);
      	m_apply = new Button("Apply");
      	m_reset = new Button("Reset");
      	m_applyRotation = new Button("Apply Rotation");
      	m_apply.addActionListener(this);
      	m_reset.addActionListener(this);
      	m_applyRotation.addActionListener(this);

      	add(m_apply);
      	add(m_reset);

      	m_degrees = new PuDouble("Degrees");
		m_degrees.setDefBounds(-180,180,0.1,1);
		m_degrees.addUpdateListener(this);
		m_degrees.init();
		add(m_degrees.getInfoPanel());

		m_axisRadioButton = new CheckboxGroup();
		m_yAxis = new Checkbox("Y-axis", m_axisRadioButton, true);
        m_xAxis = new Checkbox("X-axis", m_axisRadioButton, false);
        m_zAxis = new Checkbox("Z-axis", m_axisRadioButton, false);
        add(m_yAxis);
        add(m_xAxis);
        add(m_zAxis);

        add(m_applyRotation);

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
		else if(source == m_itClosestPointButton)
		{
			m_registration.iterativeClosestPoint(m_plane.getState());
		}
		else if(source == m_apply)
		{
			applyTransformation();
		}
		else if(source == m_reset)
		{
			m_r1c1.setText("1");
			m_r1c2.setText("0");
			m_r1c3.setText("0");
			m_r2c1.setText("0");
			m_r2c2.setText("1");
			m_r2c3.setText("0");
			m_r3c1.setText("0");
			m_r3c2.setText("0");
			m_r3c3.setText("1");
			applyTransformation();
		}
		else if(m_applyRotation == source)
		{
			double degrees = m_degrees.getValue();
			double radians = degrees * Math.PI / 180;
			double cos = Math.cos(radians);
			double sin = Math.sin(radians);
			if(m_yAxis.getState())
			{
				m_r1c1.setText(Double.toString(cos));
				m_r1c2.setText(Double.toString(-sin));
				m_r1c3.setText("0");
				m_r2c1.setText(Double.toString(sin));
				m_r2c2.setText(Double.toString(cos));
				m_r2c3.setText("0");
				m_r3c1.setText("0");
				m_r3c2.setText("0");
				m_r3c3.setText("1");
			}
			else if(m_xAxis.getState())
			{
				m_r1c1.setText("1");
				m_r1c2.setText("0");
				m_r1c3.setText("0");
				m_r2c1.setText("0");
				m_r2c2.setText(Double.toString(cos));
				m_r2c3.setText(Double.toString(-sin));
				m_r3c1.setText("0");
				m_r3c2.setText(Double.toString(sin));
				m_r3c3.setText(Double.toString(cos));
			}
			else if(m_zAxis.getState())
			{
				m_r1c1.setText(Double.toString(cos));
				m_r1c2.setText("0");
				m_r1c3.setText(Double.toString(sin));
				m_r2c1.setText("0");
				m_r2c2.setText("1");
				m_r2c3.setText("0");
				m_r3c1.setText(Double.toString(-sin));
				m_r3c2.setText("0");
				m_r3c3.setText(Double.toString(cos));
			}
			applyTransformation();
		}
	}

	private void applyTransformation()
	{
		try
		{
			PdMatrix transformation = new PdMatrix(3, 3);
			transformation.setEntry(0, 0, Double.parseDouble(m_r1c1.getText()));
			transformation.setEntry(0, 1, Double.parseDouble(m_r1c2.getText()));
			transformation.setEntry(0, 2, Double.parseDouble(m_r1c3.getText()));
			transformation.setEntry(1, 0, Double.parseDouble(m_r2c1.getText()));
			transformation.setEntry(1, 1, Double.parseDouble(m_r2c2.getText()));
			transformation.setEntry(1, 2, Double.parseDouble(m_r2c3.getText()));
			transformation.setEntry(2, 0, Double.parseDouble(m_r3c1.getText()));
			transformation.setEntry(2, 1, Double.parseDouble(m_r3c2.getText()));
			transformation.setEntry(2, 2, Double.parseDouble(m_r3c3.getText()));
			System.out.println("Calling transform");
			m_registration.transform(m_registration.m_surfP, transformation);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	/**
	 * Get information which bottom buttons a dialog should create
	 * when showing this info panel.
	 */
	protected int getDialogButtons()		{
		return PsDialog.BUTTON_OK;
	}
}
