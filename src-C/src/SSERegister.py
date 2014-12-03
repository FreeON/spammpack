class SSERegister:
  """A class that takes care of returning an available SSE register to the
  caller. This is how we rotate through the 16 available registers on an x86_64
  CPU."""

  variables = {}
  registerPool = [ "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4", "%xmm5",
      "%xmm6", "%xmm7", "%xmm8", "%xmm9", "%xmm10", "%xmm11", "%xmm12",
      "%xmm13", "%xmm14", "%xmm15" ]

  def __init__ (self, log, name, registerName = None):
    """Get a free register for the named variable. The register needs to be
    released again to become free again."""

    self.log = log

    if name in SSERegister.variables:
      self.log.error("The variable %s is already assigned a register" % (name))
      sys.exit(1)

    if registerName and not registerName in SSERegister.registerPool:
      self.log.error("The requested register %s is not available anymore" % (registerName))
      sys.exit(1)

    if registerName:
      self.name = name
      self.register = registerName
      SSERegister.registerPool.remove(registerName)
      SSERegister.variables[name] = registerName
    else:
      if len(SSERegister.registerPool) > 0:
        self.name = name
        self.register = SSERegister.registerPool.pop(0)
        SSERegister.variables[name] = self.register
      else:
        self.log.error("no registers left")
        sys.exit(1)

    self.log.debug("assigned %s --> %s" % (self.register, self.name))
    self.log.debug("registerPool: %s" % (SSERegister.registerPool))
    self.log.debug("variables: %s" % (SSERegister.variables))

  def release (self):
    """Release a register back to the register pool."""

    if not self.name in SSERegister.variables:
      self.log.error("The variable %s is not assigned a register" % (self.name))
      sys.exit(1)

    SSERegister.registerPool.append(self.register)
    del SSERegister.variables[self.name]

    self.log.debug("released register %s assigned to variable %s" % (self.register, self.name))
    self.log.debug("registerPool: %s" % (SSERegister.registerPool))
    self.log.debug("variables: %s" % (SSERegister.variables))

  def getName (self):
    """Return the name of the variable."""
    return self.name

  def getRegister (self):
    """Return the register assigned to this variable."""
    return self.register

  def __str__ (self):
    """Use this variable in the code."""

    if not self.name in SSERegister.variables:
      self.log.error("The variable %s is not assigned a register" % (self.name))
      sys.exit(1)
    return self.register
