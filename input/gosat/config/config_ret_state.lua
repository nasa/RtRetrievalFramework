------------------------------------------------------------
--- Loads state apriori from a previous L2 retrieval output
--- file.
------------------------------------------------------------

require "gosat_base_config"

config = GosatBaseConfig:new()

-- Turn on diagnostics for this specicial type of config
config.diagnostic = true

retrieved_initial_guess = ConfigCommon.oco_forward_model:new()
function retrieved_initial_guess:initial_guess()
   -- Load retrieval results from L2 output file
   l2_retrieval_file = os.getenv("l2_retrieval_file")

   if not l2_retrieval_file then
      error('The "l2_retrieval_file" environmental variable must be defined')
   end

   l2_ret = HdfFile(l2_retrieval_file)
   ret_sv = l2_ret:read_double_2d("/RetrievedStateVector/state_vector_result")(0, Range.all())

   fm_ig = ConfigCommon.oco_forward_model.initial_guess(self)

   -- After state vector structure has been built update with value from retrieved file
   if ret_sv:rows() ~= fm_ig:initial_guess():rows() then
      error("Retrieved file's state vector length: " .. ret_sv:rows() .. " does not match the configured state vector length: " .. fm_ig:initial_guess():rows())
   end

   -- Create a new initial guess which copies from the retrieved
   -- state vector and the set up covariance matrix
   ret_ig = InitialGuessValue()
   ret_ig.apriori = ret_sv
   ret_ig.initial_guess = ret_sv
   ret_ig.apriori_covariance = fm_ig:apriori_covariance()

   -- Reset the FM initial guess object
   fm_ig = CompositeInitialGuess()
   fm_ig:add_builder(ret_ig)

   return fm_ig
end
config.fm.creator = retrieved_initial_guess

-- Apply configuration
config:do_config()
