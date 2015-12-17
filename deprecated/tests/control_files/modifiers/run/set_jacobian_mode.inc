begin FILE
  source_filename = oco_l2.run

  begin PICK
      section  = CONTROL
      keyword  = jacobian_mode
      template = <use_jacobian_mode>
  end PICK

end FILE
