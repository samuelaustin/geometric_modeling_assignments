package workshop;

import java.awt.Button;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.Label;
import java.awt.CheckboxGroup;
import java.awt.Checkbox;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import jv.number.PuDouble;
import jv.object.PsConfig;
import jv.object.PsDialog;
import jv.object.PsUpdateIf;
import jvx.project.PjWorkshop_IP;

public class Smoothing_IP extends PjWorkshop_IP implements ActionListener {

	protected 	Button 			m_applySmoothing;
	protected	CheckboxGroup 	m_methodRadioButton;
	protected	Checkbox 		m_average;
	protected	Checkbox 		m_explicit;
	protected	Checkbox 		m_implicit;
	protected 	PuDouble 		m_timeStep;
	
	Smoothing m_smooth;
	
	public Smoothing_IP() {
		super();
		if(getClass() == Smoothing_IP.class)
			init();
	}
	
	public void init() {
		super.init();
		setTitle("Smoothing tool");
	}
	
	public String getNotice() {
		return "Tool for smoothing triangle meshes.";
	}
	
	public void setParent(PsUpdateIf parent) {
		super.setParent(parent);
		m_smooth = (Smoothing)parent;

		m_methodRadioButton = new CheckboxGroup();
		m_average = new Checkbox("Average", m_methodRadioButton, true);
        m_explicit = new Checkbox("Explicit", m_methodRadioButton, false);
        m_implicit = new Checkbox("Implicit", m_methodRadioButton, false);
        add(m_average);
        add(m_explicit);
        add(m_implicit);

        m_timeStep = new PuDouble("Timesteps");
		m_timeStep.setDefBounds(0,10,0.001,1);
		m_timeStep.addUpdateListener(this);
		m_timeStep.init();
		add(m_timeStep.getInfoPanel());
		
		m_applySmoothing = new Button("Apply smoothing");
		m_applySmoothing.addActionListener(this);
		add(m_applySmoothing);

		validate();
	}
	
	
	public boolean update(Object event) {
		return super.update(event);
	}
	
	/**
	 * Handle action events fired by buttons etc.
	 */
	public void actionPerformed(ActionEvent event) {
		Object source = event.getSource();
		if (source == m_applySmoothing) {
			if(m_average.getState())
			{
				m_smooth.iteratedAveraging(m_timeStep.getValue());
			}
			else if(m_explicit.getState())
			{
				m_smooth.explicitMeanCurvatureFlow(m_timeStep.getValue());
			}
			else if(m_implicit.getState())
			{
				m_smooth.implicitMeanCurvatureFlow(m_timeStep.getValue());
			}
			return;
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